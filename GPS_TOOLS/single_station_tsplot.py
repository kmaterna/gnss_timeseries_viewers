# Makes a basic python plot of time series position
# Takes in a namedTuple collection with the data. 
# Makes a basic plot. 

import numpy as np
import matplotlib.pyplot as plt 
import collections
import sys, os
import datetime as dt 
from scipy import signal
import gps_io_functions
import gps_ts_functions

# For reference of how this gets returned from the read functions.
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Parameters = collections.namedtuple("Parameters",['station','pbo_filename','unr_filename','unr_coords','datasource','outliers_remove', 'outliers_def',
	'earthquakes_remove','earthquakes_dir','offsets_remove','offsets_dir','reference_frame','seasonals_remove', 'seasonals_type','fit_table','grace_dir']);

# Types of seasonal options: 
#    fit: fits seasonals and linear trend by least squares inversion.
#   noel: uses noel's fits of inter-SSE velocities and seasonal terms.
#  notch: removes the 1-year and 6-month components by notch filter.
#  grace: uses GRACE loading model interpolated between monthly points where available, and linear inversion where not available.
#    stl: not supported yet. 


def view_single_station(station_name, offsets_remove=1, earthquakes_remove=0, outliers_remove=0, seasonals_remove=0, seasonals_type='fit',datasource='pbo'):
	MyParams=configure(station_name, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, seasonals_type, datasource);
	myData = input_data(MyParams.station, MyParams.pbo_filename, MyParams.unr_filename, MyParams.unr_coords, MyParams.datasource);
	print(myData.coords);
	[updatedData, detrended]=compute(myData,MyParams);
	single_ts_plot(updatedData,detrended,MyParams);


# Mid stream: Either have Configure determine which input file format, or have input_data return the datasource. 
# It turns out that the UNR offsets have a different format, which doesn't give you the value of the offset, but instead only gives you the time. 
# This is gonna get a little annoying. 


# -------------- CONFIGURE ------------ # 
def configure(station, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, seasonals_type, datasource):
	pbo_filename="../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos"
	unr_filename="../GPS_POS_DATA/UNR_Data/"+station+".NA12.tenv3"
	unr_coords="../GPS_POS_DATA/UNR_DATA/UNR_coords_july2018.txt"
	datasource=determine_datasource(datasource, pbo_filename, unr_filename);  # tell us which directory to use. 
	earthquakes_dir="../GPS_POS_DATA/PBO_Event_Files/"
	offsets_dir="../GPS_POS_DATA/Offsets/"
	fit_table="../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt"
	grace_dir="../GPS_POS_DATA/GRACE_loading_model/"
	outliers_def       = 15.0;  # mm away from average. 
	reference_frame    = 0;
	MyParams=Parameters(station=station, pbo_filename=pbo_filename, unr_filename=unr_filename, unr_coords=unr_coords, datasource=datasource, outliers_remove=outliers_remove, 
		outliers_def=outliers_def, earthquakes_remove=earthquakes_remove, earthquakes_dir=earthquakes_dir, offsets_remove=offsets_remove, offsets_dir=offsets_dir, 
		reference_frame=reference_frame, seasonals_remove=seasonals_remove, seasonals_type=seasonals_type, fit_table=fit_table,grace_dir=grace_dir);
	print("------- %s --------" %(station));
	print("Viewing station %s, earthquakes_remove=%d, outliers_remove=%d, seasonals_remove=%d" % (station, earthquakes_remove, outliers_remove, seasonals_remove) );
	return MyParams;

def determine_datasource(input_datasource, pbo_filename, unr_filename):
	if input_datasource=='pbo' and os.path.isfile(pbo_filename):
		print("Using PBO file as input data. ");
		datasource='pbo'
	elif input_datasource=='pbo' and os.path.isfile(unr_filename):
		print("Using UNR as input data because PBO data file doesn't exist");
		datasource='unr';
	elif input_datasource=='unr' and not os.path.isfile(unr_filename):
		print("Error! Cannot find file in UNR database.");
		sys.exit(1);
	return datasource;


# ----------- INPUTS ---------------- # 

def input_data(station_name, pbo_filename, unr_filename, unr_coords, datasource):
	if datasource=='pbo':
		[myData]=gps_io_functions.read_pbo_pos_file(pbo_filename);  # PBO data format
	if datasource=='unr':
		[myData]=gps_io_functions.read_UNR_magnet_file(unr_filename, unr_coords);  # UNR data format
	return myData;


# -------------- COMPUTE ------------ # 
def compute(myData, MyParams):
	newData=myData; 
	if MyParams.offsets_remove==1:  # First step: remove offsets and earthquakes
		newData=gps_ts_functions.remove_offsets(newData, MyParams.offsets_dir, MyParams.datasource);
	if MyParams.outliers_remove==1:  # Second step: remove outliers
		newData=gps_ts_functions.remove_outliers(newData, MyParams.outliers_def);
	if MyParams.earthquakes_remove==1:
		newData=gps_ts_functions.remove_earthquakes(newData, MyParams.earthquakes_dir, MyParams.datasource);

	trend_out=gps_ts_functions.make_detrended_option(newData, MyParams.seasonals_remove, MyParams.seasonals_type, MyParams);
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
	if MyParams.seasonals_type=="noel":
		savename=savename+"_noelfits"
		title_name=title_name+' by interSSE data'
	if MyParams.seasonals_type=="notch":
		savename=savename+"_notch"
		title_name=title_name+' by notch filter'
	if MyParams.seasonals_type=="grace":
		savename=savename+"_grace"
		title_name=title_name+' by GRACE model'		
	savename=savename+"_ts.jpg"

	axarr[0].set_title(title_name);

	plt.savefig(savename,dpi=dpival);
	return;

