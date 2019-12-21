# The purpose of this script is to compute a common mode term (average of the time series)
# We save that as a time series format
# We plot the vertical time series with and without the common mode term. 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import scipy.ndimage
import collections
import subprocess
import sys
import datetime as dt 
import gps_io_functions
import gps_input_pipeline
import gps_ts_functions
import gps_seasonal_removals
import stations_within_radius
import offsets
import remove_ets_events
import pygmt


Parameters=collections.namedtuple("Parameters",['expname','proc_center','center','radius','stations','distances','blacklist','outdir', 'outname']);


def driver():
	myparams = configure();
	[dataobj_list, offsetobj_list, eqobj_list, paired_distances] = gps_input_pipeline.multi_station_inputs(myparams.stations, myparams.blacklist, myparams.proc_center, myparams.distances);
	[common_mode, raw_objects, cmr_objects] = compute(dataobj_list, offsetobj_list, eqobj_list);
	# vertical_filtered_plots(cmr_objects, paired_distances, myparams);
	# pygmt_map(stage2_objects,myparams);
	return;


def configure():
	# center=[-125.134, 40.829]; expname='Mend'; radius = 120; # Keeper for the 2014 earthquake
	# center=[-125.134, 40.829]; expname='Mend'; radius = 90; # 
	# center=[-124, 40.5]; expname='Humboldt'; radius = 60; # 
	# center=[-122.5, 40.5]; expname='Chico'; radius = 75; # 
	# center=[-124.0, 38.0];     expname='Nbay'; radius = 125; 
	# center=[-119.0, 34.5];     expname='SoCal';  radius = 25; # km
	# center=[-116.0, 34.5];     expname='Mojave';  radius = 35; # km
	# center=[-117.5, 35.5];     expname='ECSZ';  radius = 50; # km
	# center=[-119.0, 37.7];     expname='LVC';  radius = 30; # km
	# center=[-115.5, 32.85]; expname='SSGF'; radius = 20; 
	center=[-115.5, 33]; expname='SSGF'; radius = 15; 

	proc_center='cwu';   # WHICH DATASTREAM DO YOU WANT?

	stations, distances = stations_within_radius.get_stations_within_radius(center, radius, network=proc_center);
	blacklist=["P316","P170","P158","TRND","P203","BBDM","KBRC","RYAN","BEAT","CAEC","MEXI","BOMG"];  # This is global, just keeps growing
	outdir=expname+"_"+proc_center
	subprocess.call(["mkdir","-p",outdir],shell=False);
	outname=expname+"_"+str(center[0])+"_"+str(center[1])+"_"+str(radius)
	myparams=Parameters(expname=expname, proc_center=proc_center, center=center, radius=radius, stations=stations, distances=distances, blacklist=blacklist, outdir=outdir, outname=outname);
	return myparams;

def compute(dataobj_list, offsetobj_list, eqobj_list):

	detrended_objects=[]; raw_objects=[]; cmr_objects = [];  # common-mode-removed objects
	# Objects with no earthquakes or seasonals
	for i in range(len(dataobj_list)):
		# Remove the steps and earthquakes
		newobj=offsets.remove_offsets(dataobj_list[i], offsetobj_list[i]);
		newobj=offsets.remove_offsets(newobj,eqobj_list[i]);
		newobj=gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq');
		detrended_objects.append(newobj);

	common_mode_obj = define_common_mode(detrended_objects);

	return [common_mode_obj, raw_objects, cmr_objects];


def define_common_mode(detrended_objects):
	common_mode_obj=[];
	total_dtarray, data_array_dE, data_array_dN, data_array_dU = make_ND_arrays(detrended_objects)

	plt.figure()
	# for i in range(len(detrended_objects)):
	# 	plt.plot(detrended_objects[i].dtarray,detrended_objects[i].dU);
	# 	break;
	for i in range(len(detrended_objects)):
		plt.plot(total_dtarray, data_array_dE[:,i],linewidth=0.5);
	plt.savefig("temp.png");

	# YAY WE REARRANGED!  NOW WE JUST HAVE TO TAKE A MEAN/MEDIAN AND MAKE THE COMMON MODE! 
	return common_mode_obj;



def make_ND_arrays(object_list):
	# Make a 2D array for each of the components. 
	first_date = object_list[0].dtarray[0];
	last_date = object_list[0].dtarray[-1];
	total_dtarray=[];
	# Get the list of 
	for i in range(len(object_list)):
		if object_list[i].dtarray[0]<first_date:
			first_date=object_list[i].dtarray[0];
		if object_list[i].dtarray[-1]>last_date:
			last_date=object_list[i].dtarray[-1];
	total_dtarray.append(first_date);

	for i in range(1,10000):
		if total_dtarray[-1]<last_date:
			total_dtarray.append(first_date+dt.timedelta(days=i));

	data_array_dE = np.full([len(total_dtarray), len(object_list)], np.nan);
	data_array_dN = np.full([len(total_dtarray), len(object_list)], np.nan);
	data_array_dU = np.full([len(total_dtarray), len(object_list)], np.nan);

	# Here's a bit slow. 
	for i in range(len(total_dtarray)):  # for each day... 
		for j in range(len(object_list)):  # check each station
			indices = [d == total_dtarray[i] for d in object_list[j].dtarray ];   # select which day you're doing
			# A boolean list. 
			if True in indices:
				data_array_dE[i][j] = object_list[j].dE[indices]  # please don't have a duplicated day...
				data_array_dN[i][j] = object_list[j].dN[indices]  # please don't have a duplicated day...
				data_array_dU[i][j] = object_list[j].dU[indices] # please don't have a duplicated day...

	return total_dtarray, data_array_dE, data_array_dN, data_array_dU;


def vertical_filtered_plots(dataobj_list, distances, myparams):

	plt.figure(figsize=(15,8),dpi=160);
	label_date=dt.datetime.strptime("20200215","%Y%m%d");
	start_time_plot=dt.datetime.strptime("20050101","%Y%m%d");
	end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d");

	spacing=40;
	EQtimes, labeltimes, labels, closest_station, farthest_station=configure_beautiful_plots(myparams.expname, distances);
	color_boundary_object=matplotlib.colors.Normalize(vmin=closest_station,vmax=farthest_station, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# Vertical
	for i in range(len(dataobj_list)):
		umean=np.mean(dataobj_list[i].dU);  # start at the mean. 
		# umean=dataobj_list[i].dU[0];  # start at the beginning
		line_color=custom_cmap.to_rgba(distances[i]);
		# l1 = plt.gca().plot(dataobj_list[i].dtarray,dataobj_list[i].dU-umean,linestyle='solid',linewidth=0,marker='.',color='red' ); # for debugging the filter
		udata=scipy.ndimage.median_filter(dataobj_list[i].dU-umean,size=365);
		l1 = plt.gca().plot(dataobj_list[i].dtarray,udata,linestyle='solid',linewidth=1,color=line_color );
		# plt.gca().text(label_date,offset,dataobj_list[i].name,fontsize=9,color=line_color);
	plt.gca().set_xlim(start_time_plot,end_time_plot);
	bottom,top=plt.gca().get_ylim();
	for i in range(len(EQtimes)):
		plt.gca().plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	plt.gca().set_ylabel("Filtered Vertical (mm)");
	plt.gca().set_title("Filtered Vertical GPS Time Series")
	plt.gca().grid(True)

	custom_cmap.set_array(range(int(closest_station),int(farthest_station)));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilometers from center');

	# new axis for plotting the map of california
	ax=plt.axes(position=[0.8,0.1,0.1,0.2],xticklabels=[],yticklabels=[]);
	[ca_lons,ca_lats]=np.loadtxt('../california_bdr',unpack=True);
	ax.plot(ca_lons,ca_lats,'k');
	for i in range(len(dataobj_list)):
		ax.plot(dataobj_list[i].coords[0],dataobj_list[i].coords[1],'.g',markersize=0.6);

	# new axis for extra labels
	ax=plt.axes(position=[0.8,0.4,0.1,0.2],xticklabels=[],yticklabels=[]);
	ax.set_ylim([-1,1])
	ax.set_xlim([-0.1,1])
	ax.text(0,0.75,myparams.proc_center+" data");
	ax.text(0,0.37,myparams.center);
	ax.text(0,0,str(myparams.radius)+" km radius");

	plt.savefig(myparams.outdir+"/"+myparams.outname+'_TS_'+"vertical_filt"+'.jpg');
	plt.close();
	print("Vertical plot created.");
	return;



