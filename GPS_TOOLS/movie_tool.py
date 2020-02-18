# Python viewing to see many stations' deviation from normal trend in movie form
# Step 1: Make proper directories
# Step 2: Make a movie in GMT. 

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.ndimage
import collections
import subprocess, sys
import datetime as dt 
import outputs_gps_stacks


Parameters=collections.namedtuple("Parameters",['expname','proc_center','refframe','center','radius','stations','distances','blacklist','outdir', 'outname']);
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


def movie_driver(ts_objects, distances, myparams):
	myparams = configure_movie(myparams);
	movie_ts_objects = turn_into_movie_ts(ts_objects);
	write_outputs(movie_ts_objects, myparams);
	outputs_gps_stacks.horizontal_full_ts(movie_ts_objects, distances, myparams, "");
	outputs_gps_stacks.vertical_full_ts(movie_ts_objects, distances, myparams, "");
	print("Movie data printed");
	return;

def configure_movie(myparams):
	outdir_lower=myparams.outdir+"/"+myparams.expname+"_"+str(myparams.center[0])+"_"+str(myparams.center[1])+"_"+str(myparams.radius)
	subprocess.call(["mkdir","-p",outdir_lower],shell=False);
	myparams=Parameters(expname=myparams.expname, proc_center=myparams.proc_center, refframe=myparams.refframe, center=myparams.center, radius=myparams.radius, stations=myparams.stations, distances=myparams.distances, blacklist=myparams.blacklist, outdir=outdir_lower, outname=myparams.outname);
	return myparams;

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
	# start_time = dt.datetime.strptime("2006-01-05","%Y-%m-%d");
	# end_time = dt.datetime.strptime("2019-11-15","%Y-%m-%d");
	start_time = dt.datetime.strptime("2010-04-05","%Y-%m-%d");
	end_time = dt.datetime.strptime("2015-03-15","%Y-%m-%d");	
	number_of_samples = int((end_time-start_time).days / day_interval) + 1;
	dt_downsample = [start_time+dt.timedelta(days=n*day_interval) for n in range(number_of_samples)];	
	return dt_downsample; 

def interval_downsample(dtarray, dE, dN, dU, dt_downsample):
	# This is the key math of the program. 
	dt_downsample_keep=[]; dE_downsample=[]; dN_downsample=[]; dU_downsample=[]; Se_downsample=[]; Sn_downsample=[]; Su_downsample=[];

	# Filter if desired
	edata=scipy.ndimage.median_filter(dE,size=60);
	ndata=scipy.ndimage.median_filter(dN,size=60);
	udata=scipy.ndimage.median_filter(dU,size=60);
	
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

def write_outputs(objects, myparams):

	outfile0 = myparams.outdir+"/temp.txt";
	outfile1 = myparams.outdir+"/data_out.txt";
	outfile2 = myparams.outdir+"/dates_out.txt";

	unique_dates=[];
	ofile=open(outfile0,'w');  # write the displacement at each station at each day. 
	for i in range(len(objects)):
		for j in range(len(objects[i].dtarray)):
			datestring = dt.datetime.strftime(objects[i].dtarray[j], "%Y-%m-%d");
			if datestring not in unique_dates:
				unique_dates.append(datestring);
			ofile.write("%s %f %f %s %f %f %f \n" % (objects[i].name, objects[i].coords[0], objects[i].coords[1], datestring, objects[i].dE[j], objects[i].dN[j], objects[i].dU[j]) );
	ofile.close();

	# Make an average displacement value for each day
	# format: station, lon, lat, date, e, n, u, e_avg, n_avg, u_avg
	ofile=open(outfile2,'w');
	edata_means=[]; ndata_means=[]; udata_means=[];
	for item in unique_dates:
		table = subprocess.check_output("grep "+item+" "+outfile0,shell=True);
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
	ifile=open(outfile0,'r');
	for line in ifile:
		date_item = line.split()[3];
		myindex=unique_dates.index(date_item);
		ofile.write(line[0:-1])
		ofile.write("%f %f %f \n" % (edata_means[myindex], ndata_means[myindex], udata_means[myindex]) );
	ofile.close();

	return;

