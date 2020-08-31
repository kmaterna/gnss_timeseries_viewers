

import numpy as np
import matplotlib.pyplot as plt 
import datetime as dt 
import gps_io_functions
import gps_input_pipeline
import gps_ts_functions
import gps_seasonal_removals
import stations_within_radius
import offsets

def driver():
	[stations, distances, EQtimes, proc_center, expname] = configure();
	[dataobj_list, offsetobj_list, eqobj_list] = inputs(stations, proc_center);
	[detrended_objects] = compute(dataobj_list, offsetobj_list, eqobj_list, distances, EQtimes);	
	plots(detrended_objects, EQtimes, expname, proc_center);
	# vertical_plots(stage2_objects, sorted_distances, EQtimes, expname, "vertical");
	return;

def configure():
	# EQcoords=[-125.134, 40.829]; expname='Mend'; radius = 160; 
	# EQcoords=[-124.0, 38.0];     expname='Nbay'; radius = 125; 
	# EQcoords=[-119.0, 37.7];     expname='LVC';  radius = 40; # km
	EQcoords=[-118.767, 36.231];     expname='Sierra';  radius = 70; # km
	EQtimes = [];  # What black lines do you want added to the figure? 
	EQtimes.append(dt.datetime.strptime("20140310", "%Y%m%d"));  # starts with the most important one
	EQtimes.append(dt.datetime.strptime("20050615", "%Y%m%d"));  # other earthquakes added to the figure
	EQtimes.append(dt.datetime.strptime("20100110", "%Y%m%d"));
	EQtimes.append(dt.datetime.strptime("20161208", "%Y%m%d"));
	
	proc_center='gldas';   # WHICH DATASTREAM DO YOU WANT?
	stations, lons, lats, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius);
	return [stations, distances, EQtimes, proc_center, expname];


def inputs(station_names, proc_center):  # Returns a list of objects for time series data, offsets, and earthquakes
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[];
	for station_name in station_names:
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, proc_center,, "NA";
		if myData==[]:
			continue;
		if myData.dtarray[-1]>dt.datetime.strptime("20140310","%Y%m%d") and myData.dtarray[0]<dt.datetime.strptime("20100310","%Y%m%d"):  
		# kicking out the stations that end early or start late. 
			dataobj_list.append(myData);
			offsetobj_list.append(offset_obj);
			eqobj_list.append(eq_obj);
	return [dataobj_list, offsetobj_list, eqobj_list];


def compute(dataobj_list, offsetobj_list, eqobj_list, distances, EQtimes):
	# Detrended objects
	detrended_objects=[];
	for i in range(len(dataobj_list)):
		newobj=offsets.remove_offsets(dataobj_list[i], offsetobj_list[i]);
		newobj=offsets.remove_offsets(newobj,eqobj_list[i]);
		newobj=gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq');
		detrended_objects.append(newobj);
	return [detrended_objects];


def plots(detrended_objects, EQtimes, expname, exptype):

	f1 = plt.figure();
	for i in range(len(detrended_objects)):
		plt.gca().plot_date(detrended_objects[i].dtarray, detrended_objects[i].dE,'--',color='gray');
	plt.gca().grid(True);
	for i in range(len(EQtimes)):
		plt.gca().plot_date([EQtimes[i], EQtimes[i]],[-5,5],'--k');
	plt.gca().set_ylim([-5, 5]);
	starttime=dt.datetime.strptime("20060101","%Y%m%d");
	endtime = dt.datetime.strptime("20190101","%Y%m%d");
	plt.gca().set_xlim([starttime, endtime]);
	plt.gca().set_ylabel('East (mm)');

	# # new axis for plotting the map of california
	ax=plt.axes([0.15,0.70,0.2,0.2],xticklabels=[],yticklabels=[]);
	[ca_lons,ca_lats]=np.loadtxt('california_bdr',unpack=True);
	ax.plot(ca_lons,ca_lats,'k');
	for i in range(len(detrended_objects)):
		ax.plot(detrended_objects[i].coords[0],detrended_objects[i].coords[1],'.g',markersize=0.6);

	plt.savefig('Outputs/'+exptype+'_'+expname+'_East.png');


	f1 = plt.figure();
	for i in range(len(detrended_objects)):
		plt.gca().plot_date(detrended_objects[i].dtarray, detrended_objects[i].dU,'--',color='gray');
	plt.gca().grid(True);
	for i in range(len(EQtimes)):
		plt.gca().plot_date([EQtimes[i], EQtimes[i]],[-5,5],'--k');
	plt.gca().set_ylim([-12, 12]);
	starttime=dt.datetime.strptime("20060101","%Y%m%d");
	endtime = dt.datetime.strptime("20190101","%Y%m%d");
	plt.gca().set_xlim([starttime, endtime]);
	plt.gca().set_ylabel('Vertical (mm)');

	# # new axis for plotting the map of california
	ax=plt.axes([0.15,0.70,0.2,0.2],xticklabels=[],yticklabels=[]);
	[ca_lons,ca_lats]=np.loadtxt('california_bdr',unpack=True);
	ax.plot(ca_lons,ca_lats,'k');
	for i in range(len(detrended_objects)):
		ax.plot(detrended_objects[i].coords[0],detrended_objects[i].coords[1],'.g',markersize=0.6);

	plt.savefig('Outputs/'+exptype+'_'+expname+'_Up.png');
	return;


if __name__=="__main__":
	driver();