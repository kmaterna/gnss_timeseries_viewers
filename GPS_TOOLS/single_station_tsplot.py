# Makes a basic python plot of time series position
# Takes in a namedTuple collection with the data. 
# Makes a basic plot. 

import numpy as np
import matplotlib.pyplot as plt 
import collections
import sys, os
import datetime as dt 
import gps_io_functions
import gps_ts_functions
import gps_seasonal_removals
import offsets
import gps_input_pipeline

# For reference of how this gets returned from the read functions.
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Parameters = collections.namedtuple("Parameters",['station','outliers_remove', 'outliers_def',
	'earthquakes_remove','offsets_remove','reference_frame','seasonals_remove', 'seasonals_type','fit_table','grace_dir']);
Offsets    = collections.namedtuple("Offsets",['e_offsets', 'n_offsets', 'u_offsets', 'dtevts']);

# Types of seasonal options: 
#   lssq: fits seasonals and linear trend by least squares inversion.
#   noel: uses noel's fits of inter-SSE velocities and seasonal terms.
#  notch: removes the 1-year and 6-month components by notch filter.
#  grace: uses GRACE loading model interpolated between monthly points where available, and linear inversion where not available.
#    stl: only works when Matlab is open on your computer. 


def view_single_station(station_name, offsets_remove=1, earthquakes_remove=0, outliers_remove=0, seasonals_remove=0, seasonals_type='fit',datasource='pbo'):
	MyParams = configure(station_name, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, seasonals_type);
	[myData, offset_obj, eq_obj] = input_data(MyParams.station, datasource);
	[updatedData, detrended] = compute(myData, offset_obj, eq_obj, MyParams);
	single_ts_plot(updatedData,detrended,MyParams);


# -------------- CONFIGURE ------------ # 
def configure(station, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, seasonals_type):
	fit_table="../../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt"
	grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"
	outliers_def       = 15.0;  # mm away from average. 
	reference_frame    = 0;
	MyParams=Parameters(station=station, outliers_remove=outliers_remove, outliers_def=outliers_def, earthquakes_remove=earthquakes_remove, 
		offsets_remove=offsets_remove, reference_frame=reference_frame, seasonals_remove=seasonals_remove, seasonals_type=seasonals_type, 
		fit_table=fit_table, grace_dir=grace_dir);
	print("------- %s --------" %(station));
	print("Viewing station %s, earthquakes_remove=%d, outliers_remove=%d, seasonals_remove=%d" % (station, earthquakes_remove, outliers_remove, seasonals_remove) );
	return MyParams;


# ----------- INPUTS ---------------- # 
def input_data(station_name, datasource):
	[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, datasource);
	return [myData, offset_obj, eq_obj];


# -------------- COMPUTE ------------ # 
def compute(myData, offset_obj, eq_obj, MyParams):
	newData=myData; 
	if MyParams.offsets_remove==1:  # First step: remove offsets and earthquakes
		newData=offsets.remove_antenna_offsets(newData, offset_obj);
	if MyParams.outliers_remove==1:  # Second step: remove outliers
		newData=gps_ts_functions.remove_outliers(newData, MyParams.outliers_def);
	if MyParams.earthquakes_remove==1:
		newData=offsets.remove_earthquakes(newData, eq_obj);
	
	trend_out=gps_seasonal_removals.make_detrended_ts(newData, MyParams.seasonals_remove, MyParams.seasonals_type, MyParams.fit_table, MyParams.grace_dir);
	return [newData, trend_out];


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, detrended, MyParams):
	# The major figure
	dpival=500;
	plt.figure(figsize=(15,15),dpi=dpival);
	[f,axarr]=plt.subplots(3,1,sharex=True);
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
	title_name=ts_obj.name;
	if MyParams.earthquakes_remove:
		savename=savename+"_noeq";
		title_name=title_name+', without earthquakes'
	if MyParams.seasonals_remove:
		savename=savename+"_noseasons";
		title_name=title_name+', without seasonals'
	
	if MyParams.seasonals_type=="lssq":
		savename=savename+"_lssq"
		title_name=title_name+' by least squares'
	if MyParams.seasonals_type=="noel":
		savename=savename+"_noelfits"
		title_name=title_name+' by interSSE data'
	if MyParams.seasonals_type=="notch":
		savename=savename+"_notch"
		title_name=title_name+' by notch filter'
	if MyParams.seasonals_type=="grace":
		savename=savename+"_grace"
		title_name=title_name+' by GRACE model'
	if MyParams.seasonals_type=="stl":
		savename=savename+"_stl"
		title_name=title_name+' by STL'			
	savename=savename+"_ts.jpg"

	axarr[0].set_title(title_name);

	plt.savefig(savename,dpi=dpival);
	return;

