# Makes a basic python plot of time series position
# Takes in a namedTuple collection with the data. 
# Makes a basic plot. 

import numpy as np
import matplotlib.pyplot as plt 
import collections
import datetime as dt 
from scipy import signal
import gps_io_functions
import gps_ts_functions

# For reference of how this gets returned from the read functions.
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Parameters = collections.namedtuple("Parameters",['station','filename','outliers_remove', 'outliers_def','earthquakes_remove','earthquakes_dir','offsets_remove','offsets_dir','reference_frame']);


def view_single_station(station_name, offsets_remove=1, earthquakes_remove=0, outliers_remove=0):
	MyParams=configure(station_name, offsets_remove, earthquakes_remove, outliers_remove);
	[myData]=gps_io_functions.read_pbo_pos_file(MyParams.filename);
	[updatedData, detrended]=compute(myData,MyParams);
	single_ts_plot(updatedData,detrended,MyParams);


# -------------- CONFIGURE ------------ # 
def configure(station, offsets_remove, earthquakes_remove, outliers_remove):
	filename="../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos"
	earthquakes_dir="../GPS_POS_DATA/Event_Files/"
	offsets_dir="../GPS_POS_DATA/Offsets/"
	outliers_def       = 15.0;  # mm away from average. 
	reference_frame    = 0;
	MyParams=Parameters(station=station,filename=filename, outliers_remove=outliers_remove, outliers_def=outliers_def, 
		earthquakes_remove=earthquakes_remove, earthquakes_dir=earthquakes_dir, offsets_remove=offsets_remove, offsets_dir=offsets_dir, reference_frame=reference_frame);
	print("Viewing station %s, earthquakes_remove=%d, outliers_remove=%d" % (station, earthquakes_remove, outliers_remove) );
	return MyParams;


# -------------- COMPUTE ------------ # 
def compute(myData, MyParams):
	newData=myData; 
	if MyParams.offsets_remove==1:  # First step: remove offsets and earthquakes
		newData=gps_ts_functions.remove_offsets(newData, MyParams.offsets_dir);
	if MyParams.earthquakes_remove==1:
		newData=gps_ts_functions.remove_earthquakes(newData, MyParams.earthquakes_dir);	
	if MyParams.outliers_remove==1:  # Second step: remove outliers
		newData=gps_ts_functions.remove_outliers(newData, MyParams.outliers_def);
	detrended=gps_ts_functions.detrend_data(newData);  # a ts object with detrended data
	return [newData, detrended];


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, detrended, MyParams):

	# The major figure
	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].set_title(ts_obj.name);
	axarr[0].plot_date(ts_obj.dtarray, ts_obj.dE);
	axarr[0].grid('on');
	axarr[0].set_ylabel('east (mm)');
	bottom,top=axarr[0].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[0].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k');
	ax1=axarr[0].twinx();
	ax1.plot_date(detrended.dtarray, detrended.dE,'r');
	ax1.set_ylabel('detrended (mm)')

	axarr[1].plot_date(ts_obj.dtarray, ts_obj.dN);
	axarr[1].grid('on');
	axarr[1].set_ylabel('north (mm)');
	bottom,top=axarr[1].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[1].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k');	
	ax2=axarr[1].twinx();
	ax2.plot_date(detrended.dtarray, detrended.dN,'r');
	ax2.set_ylabel('detrended (mm)')
	
	axarr[2].plot_date(ts_obj.dtarray, ts_obj.dU);	
	axarr[2].grid('on');
	axarr[2].set_ylabel('vertical (mm)')
	bottom,top=axarr[2].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[2].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k');	
	ax3=axarr[2].twinx();
	ax3.plot_date(detrended.dtarray, detrended.dU,'r');
	ax3.set_ylabel('detrended (mm)')
	

	if MyParams.earthquakes_remove==1:
		plt.savefig("single_plots/"+ts_obj.name+"_ts_noeq.jpg");
	else:
		plt.savefig("single_plots/"+ts_obj.name+"_ts.jpg");

	return;

