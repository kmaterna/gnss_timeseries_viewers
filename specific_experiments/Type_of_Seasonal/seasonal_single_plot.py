#   Goal: Develop tools for removing seasonals by...
#    fit: fits seasonals and linear trend by least squares inversion.
#   noel: uses noel's fits of inter-SSE velocities and seasonal terms.
#  notch: removes the 1-year and 6-month components by notch filter.
#  grace: uses GRACE loading model interpolated between monthly points where available, and linear inversion where not available.
#    stl: Cleveland et al., 1990 method. 
#  nldas: uses NLDAS loading model
#  gldas: uses GLDAS loading model
#   lsdm: uses LSDM loading model


import numpy as np
import matplotlib.pyplot as plt 
import collections, sys
import datetime as dt 
import gps_io_functions
import gps_ts_functions
import gps_seasonal_removals
import gps_input_pipeline
import offsets

# For reference of how this gets returned from the read functions.
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Parameters = collections.namedtuple("Parameters",['station','outliers_remove', 'outliers_def',
	'earthquakes_remove','offsets_remove','reference_frame','fit_table', 'grace_dir']);


def compare_single_seasonals(station_name, offsets_remove=1, earthquakes_remove=0, outliers_remove=0, datasource='pbo'):
	MyParams=configure(station_name, offsets_remove, earthquakes_remove, outliers_remove);
	[myData, offset_obj, eq_obj] = input_data(station_name, datasource);
	[updatedData, lssq_fit, notch_filt, grace_filt, stl_filt, gldas_filt, nldas_filt, lsdm_fit]=compute(myData, offset_obj, eq_obj, MyParams);
	single_ts_plot(updatedData,lssq_fit,notch_filt,grace_filt,stl_filt, gldas_filt, nldas_filt, lsdm_fit, MyParams);


# -------------- CONFIGURE ------------ # 
def configure(station, offsets_remove, earthquakes_remove, outliers_remove):
	fit_table="../../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt"
	grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"
	outliers_def       = 15.0;  # mm away from average. 
	reference_frame    = 0;
	MyParams=Parameters(station=station,outliers_remove=outliers_remove, outliers_def=outliers_def, 
		earthquakes_remove=earthquakes_remove, offsets_remove=offsets_remove, reference_frame=reference_frame, fit_table=fit_table, grace_dir=grace_dir);
	print("------- %s --------" %(station));
	print("Viewing station %s, earthquakes_remove=%d, outliers_remove=%d " % (station, earthquakes_remove, outliers_remove) );
	return MyParams;


# ----------- INPUTS ---------------- # 
def input_data(station_name, datasource):
	[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, datasource);
	return [myData, offset_obj, eq_obj];


# -------------- COMPUTE ------------ # 
def compute(myData, offset_obj, eq_obj, MyParams):
	newData=myData; 
	if MyParams.offsets_remove==1:  # First step: remove offsets and earthquakes
		newData=offsets.remove_offsets(newData, offset_obj);
	if MyParams.outliers_remove==1:  # Second step: remove outliers
		newData=gps_ts_functions.remove_outliers(newData, MyParams.outliers_def);
	if MyParams.earthquakes_remove==1:
		newData=offsets.remove_offsets(newData, eq_obj);

	# A few different types of seasonal removal. 
	notrend=gps_seasonal_removals.make_detrended_ts(newData, 0, 'lssq');
	lssq_fit=gps_seasonal_removals.make_detrended_ts(newData, 1, 'lssq');
	notch_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'notch');
	grace_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'grace');
	stl_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'stl');
	gldas_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'gldas');
	nldas_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'nldas');
	lsdm_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'lsdm');
	# stl_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'lssq');

	return [notrend, lssq_fit, notch_filt, grace_filt, stl_filt, gldas_filt, nldas_filt, lsdm_filt];


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, lssq_fit, notch_filt, grace_filt, stl_filt, gldas_filt, nldas_filt, lsdm_filt, MyParams):
	# The major figure
	dpival=100;
	offset_val=15;
	text_val=8;
	labeldate=dt.datetime.strptime("20070101", "%Y%m%d");
	EQtimes=[];
	EQtimes.append(dt.datetime.strptime("20050615","%Y%m%d"));
	EQtimes.append(dt.datetime.strptime("20100110","%Y%m%d"));
	EQtimes.append(dt.datetime.strptime("20140310","%Y%m%d"));
	EQtimes.append(dt.datetime.strptime("20161208","%Y%m%d"));

	plt.figure(figsize=(15,15),dpi=dpival);
	plt.plot_date(ts_obj.dtarray, ts_obj.dE+0*offset_val,color='blue',markeredgecolor='black',markersize=1.5);
	plt.plot_date(lssq_fit.dtarray, lssq_fit.dE+1*offset_val,marker='D',markersize=1.0,color='red');
	plt.plot_date(notch_filt.dtarray, notch_filt.dE+2*offset_val,marker='D',markersize=1.0,color='blue');
	plt.plot_date(stl_filt.dtarray, stl_filt.dE+3*offset_val,marker='D',markersize=1.0,color='magenta');
	plt.plot_date(grace_filt.dtarray, grace_filt.dE+4*offset_val,marker='D',markersize=1.0,color='green');
	plt.plot_date(nldas_filt.dtarray, nldas_filt.dE+5*offset_val,marker='D',markersize=1.0,color='cyan');
	plt.plot_date(gldas_filt.dtarray, gldas_filt.dE+6*offset_val,marker='D',markersize=1.0,color='orange');
	plt.plot_date(lsdm_filt.dtarray, lsdm_filt.dE+7*offset_val,marker='D',markersize=1.0,color='orange');
	plt.grid(linestyle='--',linewidth=0.5);	
	plt.text(labeldate,0*offset_val+text_val,'Uncorrected',fontsize=22,color='black');
	plt.text(labeldate,1*offset_val+text_val,'LsSq Fit',fontsize=22,color='red');
	plt.text(labeldate,2*offset_val+text_val,'Notch filter',fontsize=22,color='blue');
	plt.text(labeldate,3*offset_val+text_val,'STL filter',fontsize=22,color='magenta');
	plt.text(labeldate,4*offset_val+text_val,'GRACE load model',fontsize=22,color='green');	
	plt.text(labeldate,5*offset_val+text_val,'NLDAS load model',fontsize=22,color='cyan');
	plt.text(labeldate,6*offset_val+text_val,'GLDAS load model',fontsize=22,color='orange');
	plt.text(labeldate,7*offset_val+text_val,'LSDM load model',fontsize=22,color='orange');
	bottom,top=plt.gca().get_ylim();
	for i in range(len(EQtimes)):
		plt.plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k',linewidth=1.5);
	plt.ylabel('detrended east (mm)',fontsize=22)
	plt.gca().tick_params(labelsize=22);
	title_name=ts_obj.name+" Seasonal Corrections - East";
	plt.title(title_name,fontsize=24);
	savename="single_plots/seasonal_"+ts_obj.name;
	savename=savename+"_east.jpg"
	plt.savefig(savename,dpi=dpival);

	plt.figure(figsize=(15,15),dpi=dpival);
	plt.plot_date(ts_obj.dtarray, ts_obj.dN+0*offset_val,color='blue',markeredgecolor='black',markersize=1.5);
	plt.plot_date(lssq_fit.dtarray, lssq_fit.dN+1*offset_val,marker='D',markersize=1.5,color='red');
	plt.plot_date(notch_filt.dtarray, notch_filt.dN+2*offset_val,marker='D',markersize=1.5,color='blue');
	plt.plot_date(stl_filt.dtarray, stl_filt.dN+3*offset_val,marker='D',markersize=1.0,color='magenta');
	plt.plot_date(grace_filt.dtarray, grace_filt.dN+4*offset_val,marker='D',markersize=1.0,color='green');
	plt.plot_date(nldas_filt.dtarray, nldas_filt.dN+5*offset_val,marker='D',markersize=1.0,color='cyan');
	if len(gldas_filt.dN)>1:
		plt.plot_date(gldas_filt.dtarray, gldas_filt.dN+6*offset_val,marker='D',markersize=1.0,color='orange');	
	plt.grid(linestyle='--',linewidth=0.5);	
	plt.text(labeldate,0*offset_val+text_val,'Uncorrected',fontsize=22,color='black');
	plt.text(labeldate,1*offset_val+text_val,'LsSq Fit',fontsize=22,color='red');
	plt.text(labeldate,2*offset_val+text_val,'Notch filter',fontsize=22,color='blue');	
	plt.text(labeldate,3*offset_val+text_val,'STL filter',fontsize=22,color='magenta');
	plt.text(labeldate,4*offset_val+text_val,'GRACE load model',fontsize=22,color='green');	
	plt.text(labeldate,5*offset_val+text_val,'NLDAS load model',fontsize=22,color='cyan');
	plt.text(labeldate,6*offset_val+text_val,'GLDAS load model',fontsize=22,color='orange');
	plt.ylabel('detrended north (mm)',fontsize=22)
	plt.gca().tick_params(labelsize=22);
	title_name=ts_obj.name+" Seasonal Corrections - North";
	plt.title(title_name,fontsize=24);
	savename="single_plots/seasonal_"+ts_obj.name;
	savename=savename+"_north.jpg"
	plt.savefig(savename,dpi=dpival);

	plt.figure(figsize=(15,15),dpi=dpival);
	plt.plot_date(ts_obj.dtarray, ts_obj.dU+0*offset_val,color='blue',markeredgecolor='black',markersize=1.5);
	plt.plot_date(lssq_fit.dtarray, lssq_fit.dU+2*offset_val,marker='D',markersize=1.5,color='red');
	plt.plot_date(notch_filt.dtarray, notch_filt.dU+4*offset_val,marker='D',markersize=1.5,color='blue');
	plt.plot_date(stl_filt.dtarray, stl_filt.dU+6*offset_val,marker='D',markersize=1.0,color='magenta');
	plt.plot_date(grace_filt.dtarray, grace_filt.dU+8*offset_val,marker='D',markersize=1.0,color='green');	
	plt.plot_date(nldas_filt.dtarray, nldas_filt.dU+10*offset_val,marker='D',markersize=1.0,color='cyan');
	plt.plot_date(gldas_filt.dtarray, gldas_filt.dU+12*offset_val,marker='D',markersize=1.0,color='orange');
	plt.grid(linestyle='--',linewidth=0.5);	
	plt.text(labeldate,0*offset_val+text_val,'Uncorrected',fontsize=22,color='black');
	plt.text(labeldate,2*offset_val+text_val,'LsSq Fit',fontsize=22,color='red');
	plt.text(labeldate,4*offset_val+text_val,'Notch filter',fontsize=22,color='blue');	
	plt.text(labeldate,6*offset_val+text_val,'STL filter',fontsize=22,color='magenta');
	plt.text(labeldate,8*offset_val+text_val,'GRACE load model',fontsize=22,color='green');	
	plt.text(labeldate,10*offset_val+text_val,'NLDAS load model',fontsize=22,color='cyan');
	plt.text(labeldate,12*offset_val+text_val,'GLDAS load model',fontsize=22,color='orange');
	plt.ylabel('detrended vertical (mm)',fontsize=22)
	plt.gca().tick_params(labelsize=22);
	title_name=ts_obj.name+" Seasonal Corrections - Vertical";
	plt.title(title_name,fontsize=24);
	savename="single_plots/seasonal_"+ts_obj.name;
	savename=savename+"_vert.jpg"
	plt.savefig(savename,dpi=dpival);
	return;
