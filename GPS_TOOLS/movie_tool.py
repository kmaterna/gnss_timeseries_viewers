# Python viewing to see many stations' deviation from normal trend in movie form
# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and trend from the data. 
# Step 4: Make a movie in GMT. 

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.ndimage
import collections
import subprocess, sys
import datetime as dt 
import gps_io_functions
import gps_input_pipeline
import gps_ts_functions
import gps_seasonal_removals
import stations_within_radius
import outputs_gps_stacks
import offsets


Parameters=collections.namedtuple("Parameters",['expname','proc_center','refframe','center','radius','stations','distances','blacklist','outdir', 'outname']);
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

def driver():
	myparams = configure();
	[dataobj_list, offsetobj_list, eqobj_list, paired_distances] = gps_input_pipeline.multi_station_inputs(myparams.stations, myparams.blacklist, myparams.proc_center, myparams.refframe, myparams.distances);
	[detrend_objects, no_offset_objects, no_offsets_no_trends, no_offsets_no_trends_no_seasons, sorted_distances] = compute(dataobj_list, offsetobj_list, eqobj_list, paired_distances);
	movie_ts_objects = turn_into_movie_ts(no_offsets_no_trends_no_seasons);
	write_outputs(movie_ts_objects, myparams.outdir+"/data_out.txt", myparams.outdir+"/dates_out.txt");
	horizontal_full_ts(movie_ts_objects, sorted_distances, myparams, "");
	vertical_full_ts(movie_ts_objects, sorted_distances, myparams, "");
	return;

def configure():
	# center=[-115.5, 32.85]; expname='SSGF'; radius = 20; 
	# center=[-115.5, 33]; expname='SSGF'; radius = 15; 
	center=[-115.5, 33]; expname='SSGF'; radius =25; 

	proc_center='nmt';   # WHICH DATASTREAM DO YOU WANT?
	refframe = 'NA';     # WHICH REFERENCE FRAME? 

	stations, distances = stations_within_radius.get_stations_within_radius(center, radius, network=proc_center);
	blacklist=["P316","P170","P158","TRND","P203","BBDM","KBRC","RYAN","BEAT","CAEC","MEXI","BOMG","FSHB"];  # This is global, just keeps growing
	outdir_upper=expname+"_"+proc_center+"_"+refframe
	outdir_lower=outdir_upper+"/"+expname+"_"+str(center[0])+"_"+str(center[1])+"_"+str(radius)
	subprocess.call(["mkdir","-p",outdir_upper],shell=False);
	subprocess.call(["mkdir","-p",outdir_lower],shell=False);
	outname=expname+"_"+str(center[0])+"_"+str(center[1])+"_"+str(radius)
	myparams=Parameters(expname=expname, proc_center=proc_center, refframe=refframe, center=center, radius=radius, stations=stations, distances=distances, blacklist=blacklist, outdir=outdir_lower, outname=outname);
	return myparams;

def compute(dataobj_list, offsetobj_list, eqobj_list, distances):
	# This is an important function. 
	# It sorts a list of stations in latitude order, and then removes all the things you might want to remove during standard processing. 
	# Removes offsets, removes trends, removes seasonals. 
	# Returns a bunch of lists of objects, to use as you might see fit. 

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_offsets = [x for _,x in sorted(zip(latitudes_list, offsetobj_list))];  # the raw, sorted data. 
	sorted_eqs = [x for _,x in sorted(zip(latitudes_list, eqobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.
	
	detrended_objects = [];
	no_offset_objects = [];
	no_offsets_no_trends = [];
	no_offsets_no_trends_no_seasons = [];
	
	# Detrended objects (or objects with trends and no offsets; depends on what you want.)
	for i in range(len(sorted_objects)):
		newobj=gps_seasonal_removals.make_detrended_ts(sorted_objects[i], 0, 'lssq');
		detrended_objects.append(newobj);  # still has offsets, doesn't have trends
		
		newobj=offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
		newobj=offsets.remove_offsets(newobj,sorted_eqs[i]);
		no_offset_objects.append(newobj);  # still has trends, doesn't have offsets

	# Objects with no earthquakes or seasonals
	for i in range(len(dataobj_list)):

		# Remove the steps earthquakes
		newobj=offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
		newobj=offsets.remove_offsets(newobj,sorted_eqs[i]);
		newobj=gps_ts_functions.remove_outliers(newobj, 20);  # 20mm outlier definition

		# The detrended TS without earthquakes
		stage1obj=gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq');
		no_offsets_no_trends.append(stage1obj);

		# The detrended TS without earthquakes or seasonals
		stage2obj=gps_seasonal_removals.make_detrended_ts(stage1obj, 1, 'lssq');
		no_offsets_no_trends_no_seasons.append(stage2obj);

	return [detrended_objects, no_offset_objects, no_offsets_no_trends, no_offsets_no_trends_no_seasons, sorted_distances];



def turn_into_movie_ts(myobjects):
	# Here we take a raw time series, and turn it into something that we will movie-ify. 
	# For example, we might filter it, downsample it monthly, etc. 
	ds_objects=[];
	desired_downsample = get_downsample_dates();
	for i in range(len(myobjects)):
		myobject=myobjects[i];
		[dt_ds, dE_ds, dN_ds, dU_ds, Se_ds, Sn_ds, Su_ds] = interval_downsample(myobject.dtarray, myobject.dE, myobject.dN, myobject.dU, desired_downsample);
		downsampled_object = Timeseries(name=myobject.name, coords=myobject.coords, dtarray=dt_ds, dE=dE_ds, dN=dN_ds, dU=dU_ds, Se=Se_ds, Sn=Sn_ds, Su=Su_ds, EQtimes=myobject.EQtimes);
		ds_objects.append(downsampled_object);

		# f,axarr = plt.subplots(3,1, figsize=(15,15));
		# axarr[0].plot(myobject.dtarray, myobject.dE-np.nanmean(myobject.dE), 'o', markersize=1, color='black',linewidth=2);
		# axarr[0].plot(downsampled_object.dtarray, downsampled_object.dE, 'r');
		# axarr[0].set_ylabel('East (mm)');
		# axarr[0].set_title(myobject.name,fontsize=20);
		# axarr[1].plot(myobject.dtarray, myobject.dN-np.nanmean(myobject.dN), 'o', markersize=1, color='black',linewidth=2);
		# axarr[1].plot(downsampled_object.dtarray, downsampled_object.dN, 'r');
		# axarr[1].set_ylabel('North (mm)');
		# axarr[2].plot(myobject.dtarray, myobject.dU-np.nanmean(myobject.dU), 'o', markersize=1, color='black',linewidth=2);
		# axarr[2].plot(downsampled_object.dtarray, downsampled_object.dU, 'r');
		# axarr[2].set_ylabel('Up (mm)');
		# plt.savefig(downsampled_object.name+'_example_timeseries.png');

	return ds_objects; 

def get_downsample_dates():
	# Make a subsampled array of dates between start and end times. 
	day_interval = 90;   # days
	start_time = dt.datetime.strptime("2006-01-05","%Y-%m-%d");
	end_time = dt.datetime.strptime("2019-11-15","%Y-%m-%d");
	number_of_samples = int((end_time-start_time).days / day_interval) + 1;
	dt_downsample = [start_time+dt.timedelta(days=n*day_interval) for n in range(number_of_samples)];	
	return dt_downsample; 

def interval_downsample(dtarray, dE, dN, dU, dt_downsample):
	# This is the key math of the program. 
	dt_downsample_keep=[]; dE_downsample=[]; dN_downsample=[]; dU_downsample=[]; Se_downsample=[]; Sn_downsample=[]; Su_downsample=[];

	# Filter if desired
	edata=scipy.ndimage.median_filter(dE-np.nanmean(dE),size=60);
	ndata=scipy.ndimage.median_filter(dN-np.nanmean(dN),size=60);
	udata=scipy.ndimage.median_filter(dU-np.nanmean(dU),size=60);
	
	for date in dt_downsample:
		myindex=np.where(dtarray==np.datetime64(date));
		if len(myindex[0])>0:  # if the object is in the array, then we will find the exact datetime index
			myindex=myindex[0][0];
			dt_downsample_keep.append(dtarray[myindex]);
			dE_downsample.append(edata[myindex]);
			dN_downsample.append(ndata[myindex]);
			dU_downsample.append(udata[myindex]);
			Se_downsample.append(0.1);
			Sn_downsample.append(0.1);
			Su_downsample.append(0.1);
	return [dt_downsample_keep, dE_downsample, dN_downsample, dU_downsample, Se_downsample, Sn_downsample, Su_downsample]; 

def write_outputs(objects, outfile1, outfile2):
	unique_dates=[];
	ofile=open("temp.txt",'w');
	for i in range(len(objects)):
		for j in range(len(objects[i].dtarray)):
			datestring = dt.datetime.strftime(objects[i].dtarray[j], "%Y-%m-%d");
			if datestring not in unique_dates:
				unique_dates.append(datestring);
			ofile.write("%s %f %f %s %f %f %f \n" % (objects[i].name, objects[i].coords[0], objects[i].coords[1], datestring, objects[i].dE[j], objects[i].dN[j], objects[i].dU[j]) );
	ofile.close();

	# Make an average displacement value for each day
	ofile=open(outfile2,'w');
	edata_means=[]; ndata_means=[]; udata_means=[];
	for item in unique_dates:
		table = subprocess.check_output("grep "+item+" "+"temp.txt",shell=True);
		table=table.decode('utf-8');
		table_rows=table.split('\n');
		edata=[]; ndata=[]; udata=[]; 
		for line in table_rows:
			if len(line)>0:
				words = line.split();
				edata.append(float(words[4]));
				ndata.append(float(words[5]));
				udata.append(float(words[6]));
		ofile.write("%s %f %f %f" % (item, np.mean(edata), np.mean(ndata), np.mean(udata) ) );
		ofile.write("\n");
		edata_means.append(np.mean(edata))
		ndata_means.append(np.mean(ndata))
		udata_means.append(np.mean(udata))
	ofile.close();

	# Write the average displacement for each date
	ofile=open(outfile1,'w');
	ifile=open("temp.txt",'r');
	for line in ifile:
		date_item = line.split()[3];
		myindex=unique_dates.index(date_item);
		ofile.write(line[0:-1])
		ofile.write("%f %f %f \n" % (edata_means[myindex], ndata_means[myindex], udata_means[myindex]) );
	ofile.close();

	return;







