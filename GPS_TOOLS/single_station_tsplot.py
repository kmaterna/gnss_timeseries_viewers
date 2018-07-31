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
Parameters = collections.namedtuple("Parameters",['station','filename','outliers_remove', 'outliers_def',
	'earthquakes_remove','earthquakes_dir','offsets_remove','offsets_dir','reference_frame','seasonals_remove', 'fit_type','fit_table']);


def view_single_station(station_name, offsets_remove=1, earthquakes_remove=0, outliers_remove=0, seasonals_remove=0, fit_type='fit'):
	MyParams=configure(station_name, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, fit_type);
	[myData]=gps_io_functions.read_pbo_pos_file(MyParams.filename);
	print(myData.coords);
	[updatedData, detrended]=compute(myData,MyParams);
	single_ts_plot(updatedData,detrended,MyParams);


# -------------- CONFIGURE ------------ # 
def configure(station, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, fit_type):
	filename="../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos"
	earthquakes_dir="../GPS_POS_DATA/PBO_Event_Files/"
	offsets_dir="../GPS_POS_DATA/Offsets/"
	fit_table="../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt"
	outliers_def       = 15.0;  # mm away from average. 
	reference_frame    = 0;
	MyParams=Parameters(station=station,filename=filename, outliers_remove=outliers_remove, outliers_def=outliers_def, 
		earthquakes_remove=earthquakes_remove, earthquakes_dir=earthquakes_dir, offsets_remove=offsets_remove, offsets_dir=offsets_dir, 
		reference_frame=reference_frame, seasonals_remove=seasonals_remove, fit_type=fit_type, fit_table=fit_table);
	print("------- %s --------" %(station));
	print("Viewing station %s, earthquakes_remove=%d, outliers_remove=%d, seasonals_remove=%d" % (station, earthquakes_remove, outliers_remove, seasonals_remove) );
	return MyParams;


# -------------- COMPUTE ------------ # 
def compute(myData, MyParams):
	newData=myData; 
	if MyParams.offsets_remove==1:  # First step: remove offsets and earthquakes
		newData=gps_ts_functions.remove_offsets(newData, MyParams.offsets_dir);
	if MyParams.outliers_remove==1:  # Second step: remove outliers
		newData=gps_ts_functions.remove_outliers(newData, MyParams.outliers_def);
	if MyParams.earthquakes_remove==1:
		newData=gps_ts_functions.remove_earthquakes(newData, MyParams.earthquakes_dir);		
	
	if MyParams.fit_type=='fit':
		if MyParams.seasonals_remove==1:
			newData=gps_ts_functions.remove_annual_semiannual_by_fitting(newData);	
		detrended=gps_ts_functions.detrend_data_by_fitting(newData);  # a ts object with detrended data

	if MyParams.fit_type=='noel':
		if MyParams.seasonals_remove==1:
			newData=gps_ts_functions.remove_annual_semiannual_by_table(newData, MyParams.fit_table);	
		detrended=gps_ts_functions.detrend_data_by_table(newData,MyParams.fit_table);  # a ts object with detrended data
	return [newData, detrended];


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, detrended, MyParams):
	# The major figure
	dpival=500;
	plt.figure(figsize=(15,15),dpi=dpival);
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].set_title(ts_obj.name);
	axarr[0].plot_date(ts_obj.dtarray, ts_obj.dE,color='blue',markeredgecolor='black',markersize=1.5);
	axarr[0].grid(linestyle='--',linewidth=0.5);
	axarr[0].set_ylabel('east (mm)');
	bottom,top=axarr[0].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[0].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k',linewidth=1);
	ax1=axarr[0].twinx();
	ax1.plot_date(detrended.dtarray, detrended.dE,marker='D',markersize=1.0,color='red');
	ax1.set_ylabel('detrended (mm)')


	axarr[1].plot_date(ts_obj.dtarray, ts_obj.dN,color='blue',markeredgecolor='black',markersize=1.5);
	axarr[1].grid(linestyle='--',linewidth=0.5);
	axarr[1].set_ylabel('north (mm)');
	bottom,top=axarr[1].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[1].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k',linewidth=1);	
	ax2=axarr[1].twinx();
	ax2.plot_date(detrended.dtarray, detrended.dN,marker='D',markersize=1.0,color='red');
	ax2.set_ylabel('detrended (mm)')
	
	axarr[2].plot_date(ts_obj.dtarray, ts_obj.dU,color='blue',markeredgecolor='black',markersize=1.5);
	axarr[2].grid(linestyle='--',linewidth=0.5);
	axarr[2].set_ylabel('vertical (mm)')
	bottom,top=axarr[2].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[2].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k',linewidth=1);	
	ax3=axarr[2].twinx();
	ax3.plot_date(detrended.dtarray, detrended.dU,marker='D',markersize=1.0,color='red');
	ax3.set_ylabel('detrended (mm)')
	axarr[2].set_xlim([min(ts_obj.dtarray), max(ts_obj.dtarray)]);

	savename="single_plots/"+ts_obj.name;
	if MyParams.earthquakes_remove:
		savename=savename+"_noeq";
	if MyParams.seasonals_remove:
		savename=savename+"_noseasons";
	if MyParams.fit_type=="noel":
		savename=savename+"_noelfits"
	savename=savename+"_ts.jpg"

	plt.savefig(savename,dpi=dpival);
	return;

