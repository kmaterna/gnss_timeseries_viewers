# Remove data during ETS times
# Estimate offsets during those times. 

import numpy as np
import matplotlib.pyplot as plt 
import collections
import datetime as dt 
import gps_io_functions
import gps_ts_functions
import gps_seasonal_removals
import offsets
import gps_input_pipeline
import single_station_tsplot

# For reference of how this gets returned from the read functions.
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Offsets    = collections.namedtuple("Offsets",['e_offsets', 'n_offsets', 'u_offsets', 'evdts']);
Parameters = collections.namedtuple("Parameters",['station','outliers_remove', 'outliers_def',
	'earthquakes_remove','offsets_remove','reference_frame','seasonals_remove', 'seasonals_type','datasource','refframe','fit_table','grace_dir']);


def view_single_station(station_name, offsets_remove=1, earthquakes_remove=0, outliers_remove=0, seasonals_remove=0, seasonals_type='lssq',datasource='pbo',refframe='NA'):
	MyParams = configure(station_name, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, seasonals_type, datasource, refframe);
	[myData, offset_obj, eq_obj, ets_intervals] = input_data(MyParams.station, MyParams.datasource, MyParams.refframe);
	[ETS_removed, detrended1, ETS_within, detrended2] = compute(myData, offset_obj, eq_obj, MyParams, ets_intervals);
	single_ts_plot(ETS_removed,detrended1,MyParams, ets_intervals,'corrected');
	single_ts_plot(ETS_within,detrended2,MyParams, ets_intervals,'uncorrected');

# -------------- CONFIGURE ------------ # 
def configure(station, offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, seasonals_type, datasource, refframe):
	fit_table="../../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt"
	grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"
	outliers_def       = 10.0;  # mm away from average. 
	MyParams=Parameters(station=station, outliers_remove=outliers_remove, outliers_def=outliers_def, earthquakes_remove=earthquakes_remove, 
		offsets_remove=offsets_remove, reference_frame=refframe, seasonals_remove=seasonals_remove, seasonals_type=seasonals_type, 
		datasource=datasource, refframe=refframe, fit_table=fit_table, grace_dir=grace_dir);
	print("------- %s --------" %(station));
	print("Viewing station %s, earthquakes_remove=%d, outliers_remove=%d, seasonals_remove=%d, datasource=%s, refframe=%s" % (station, earthquakes_remove, outliers_remove, seasonals_remove,datasource,refframe));
	return MyParams;


# ----------- INPUTS ---------------- # 
def input_data(station_name, datasource, refframe):
	[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, datasource, refframe);
	ets_intervals=input_tremor_days();
	return [myData, offset_obj, eq_obj, ets_intervals];

def input_tremor_days():
	ifile = open("../../GPS_POS_DATA/tremor/ETS_events.txt",'r');
	ets_ints=[];
	for line in ifile:
		temp=line.split();
		t1=dt.datetime.strptime(temp[0],"%Y-%m-%d");
		t2=dt.datetime.strptime(temp[2],"%Y-%m-%d");
		ets_ints.append([t1, t2]);
	ifile.close();
	return ets_ints;

# -------------- COMPUTE ------------ # 
def compute(myData, offset_obj, eq_obj, MyParams, ets_intervals):
	# First, embed the data with the eq object (useful no matter what)
	newData=Timeseries(name=myData.name, coords=myData.coords, dtarray=myData.dtarray, dN=myData.dN, dE=myData.dE, dU=myData.dU, Sn=myData.Sn, Se=myData.Se, Su=myData.Su, EQtimes=eq_obj.evdts);
	if MyParams.offsets_remove==1:  # Remove offsets and antenna changes
		newData=offsets.remove_offsets(newData, offset_obj);
	if MyParams.outliers_remove==1:  # Remove outliers
		newData=gps_ts_functions.remove_outliers(newData, MyParams.outliers_def);
	if MyParams.earthquakes_remove==1: # Remove earthquakes
		newData=offsets.remove_offsets(newData, eq_obj);
	trend_out__uncorrected=gps_seasonal_removals.make_detrended_ts(newData, MyParams.seasonals_remove, MyParams.seasonals_type, MyParams.fit_table, MyParams.grace_dir);

	ETS_removed=remove_ETS_times(newData,ets_intervals);
	trend_out_corrected=gps_seasonal_removals.make_detrended_ts(ETS_removed, MyParams.seasonals_remove, MyParams.seasonals_type, MyParams.fit_table, MyParams.grace_dir);
	return [ETS_removed, trend_out_corrected, newData, trend_out__uncorrected];


def remove_ETS_times(ts_obj, ets_intervals):
	dtarray=[]; dE=[]; dN=[]; dU=[]; Se=[]; Sn=[]; Su=[];
	for i in range(len(ts_obj.dtarray)):
		report = True;
		for j in range(len(ets_intervals)):
			if ts_obj.dtarray[i]>=ets_intervals[j][0] and ts_obj.dtarray[i]<=ets_intervals[j][1]:
				report=False;
		if report==True and ts_obj.dtarray[i]>=dt.datetime.strptime("20120101","%Y%m%d"):
			dtarray.append(ts_obj.dtarray[i]);
			dE.append(ts_obj.dE[i]);
			dN.append(ts_obj.dN[i]);
			dU.append(ts_obj.dU[i]);
			Se.append(ts_obj.Se[i]);
			Sn.append(ts_obj.Sn[i]);
			Su.append(ts_obj.Su[i]);
	
	ts_obj_fix=Timeseries(name=ts_obj.name, coords=ts_obj.coords, dtarray=dtarray, dE=dE, dN=dN, dU=dU, Se=Se, Sn=Sn, Su=Su, EQtimes=ts_obj.EQtimes);

	# Find the offsets associated with each ETS interval
	e_offsets=[]; n_offsets=[]; u_offsets=[]; evdts=[]; 
	for i in range(len(ets_intervals)):
		e_offset=fit_offset(dtarray, dE, ets_intervals[i]);
		n_offset=fit_offset(dtarray, dN, ets_intervals[i]);
		u_offset=fit_offset(dtarray, dU, ets_intervals[i]);
		e_offsets.append(e_offset);
		n_offsets.append(n_offset);
		u_offsets.append(u_offset);
		evdts.append(ets_intervals[i][1]);
	print(e_offsets);
	offset_obj = Offsets(e_offsets=e_offsets, n_offsets=n_offsets, u_offsets=u_offsets, evdts=evdts);

	ts_obj_fix = offsets.remove_offsets(ts_obj_fix,offset_obj);
	ts_obj_new=Timeseries(name=ts_obj.name, coords=ts_obj.coords, dtarray=ts_obj_fix.dtarray, dE=ts_obj_fix.dE, dN=ts_obj_fix.dN, dU=ts_obj_fix.dU, Se=ts_obj_fix.Se, Sn=ts_obj_fix.Sn, Su=ts_obj_fix.Su, EQtimes=ts_obj_fix.EQtimes);
	return ts_obj_new;

def fit_offset(dtarray, data, interval):
	num_days=10;  
	start_dt_index=[];
	for i in range(len(dtarray)):
		deltat=dtarray[i]-interval[0];
		if abs(deltat.days)<=1:
			start_dt_index=i;
	if start_dt_index==[]:
		offset=0;
	else:
		before_mean=np.mean(data[start_dt_index-num_days:start_dt_index]);
		after_mean=np.mean(data[start_dt_index+1:start_dt_index+num_days]);
		offset=after_mean-before_mean;
	return offset;


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, detrended, MyParams, ets_intervals,name):

	label_fontsize=18;
	eq_2016 = dt.datetime.strptime("20161208","%Y%m%d");
	start_time=dt.datetime.strptime("20060101","%Y%m%d");
	end_time=dt.datetime.strptime("20180915","%Y%m%d");

	# The major figure
	dpival=500;
	[f,axarr]=plt.subplots(3,1,sharex=True,figsize=(10,7),dpi=dpival);
	axarr[0].plot_date(ts_obj.dtarray, ts_obj.dE,color='blue',markeredgecolor='black',markersize=1.5);
	axarr[0].grid(linestyle='--',linewidth=0.5);
	axarr[0].set_ylabel('east (mm)',fontsize=label_fontsize);
	bottom,top=axarr[0].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[0].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k',linewidth=0.5);
	axarr[0].plot_date([eq_2016, eq_2016], [bottom, top], '--k',linewidth=0.5);  # special plot for SSE experiment
	axarr[0].set_xlim([start_time, end_time]);
	for i in range(len(ets_intervals)):
		axarr[0].plot_date([ets_intervals[i][0],ets_intervals[i][0]],[bottom,top],'-b',linewidth=0.5);
		axarr[0].plot_date([ets_intervals[i][1],ets_intervals[i][1]],[bottom,top],'-b',linewidth=0.5);
	ax1=axarr[0].twinx();
	ax1.plot_date(detrended.dtarray, detrended.dE,marker='D',markersize=1.0,color='red');
	ax1.set_ylabel('detrended (mm)', fontsize=label_fontsize-2);
	plt.setp(axarr[0].get_xticklabels(),fontsize=label_fontsize);
	plt.setp(axarr[0].get_yticklabels(),fontsize=label_fontsize);
	plt.setp(ax1.get_yticklabels(),fontsize=label_fontsize);	


	axarr[1].plot_date(ts_obj.dtarray, ts_obj.dN,color='blue',markeredgecolor='black',markersize=1.5);
	axarr[1].grid(linestyle='--',linewidth=0.5);
	axarr[1].set_ylabel('north (mm)',fontsize=label_fontsize);
	bottom,top=axarr[1].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[1].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k',linewidth=0.5);	
	axarr[1].plot_date([eq_2016, eq_2016], [bottom, top], '--k',linewidth=0.5);  # special plot for SSE experiment
	for i in range(len(ets_intervals)):
		axarr[1].plot_date([ets_intervals[i][0],ets_intervals[i][0]],[bottom,top],'-b',linewidth=0.5);
		axarr[1].plot_date([ets_intervals[i][1],ets_intervals[i][1]],[bottom,top],'-b',linewidth=0.5);
	axarr[1].set_xlim([start_time, end_time]);		
	ax2=axarr[1].twinx();
	ax2.plot_date(detrended.dtarray, detrended.dN,marker='D',markersize=1.0,color='red');
	ax2.set_ylabel('detrended (mm)', fontsize=label_fontsize-2);
	plt.setp(axarr[1].get_xticklabels(),fontsize=label_fontsize);
	plt.setp(axarr[1].get_yticklabels(),fontsize=label_fontsize);
	plt.setp(ax2.get_yticklabels(),fontsize=label_fontsize);
	

	axarr[2].plot_date(ts_obj.dtarray, ts_obj.dU,color='blue',markeredgecolor='black',markersize=1.5);
	axarr[2].grid(linestyle='--',linewidth=0.5);
	axarr[2].set_ylabel('vertical (mm)',fontsize=label_fontsize);
	bottom,top=axarr[2].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[2].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k',linewidth=0.5);	
	ax3=axarr[2].twinx();
	ax3.plot_date(detrended.dtarray, detrended.dU,marker='D',markersize=1.0,color='red');
	ax3.set_ylabel('detrended (mm)', fontsize=label_fontsize-2);
	axarr[2].set_xlim([min(ts_obj.dtarray), max(ts_obj.dtarray)]);
	axarr[2].set_xlim([start_time, end_time]);
	plt.setp(axarr[2].get_xticklabels(),fontsize=label_fontsize);
	plt.setp(axarr[2].get_yticklabels(),fontsize=label_fontsize);
	plt.setp(ax3.get_yticklabels(),fontsize=label_fontsize);

	title=ts_obj.name+" without ETS offset";
	savename=title+"_"+name+".jpg";
	axarr[0].set_title(title,fontsize=label_fontsize+2);
	plt.savefig(savename,dpi=dpival);
	return;

if __name__=="__main__":
	station="P159"
	view_single_station(station, 
		offsets_remove=1, earthquakes_remove=1, 
		outliers_remove=1, seasonals_remove=1,seasonals_type='lssq', datasource='pbo', refframe='NA');



