#   Goal: Compare 

import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import gps_input_pipeline
import gps_io_functions
import gps_ts_functions
import offsets
import gps_seasonal_removals

def compare_single_station(station):
	Time_periods, EQs=configure_times();
	[gldas, nldas, grace, lsdm, gps] = inputs(station);
	gldas_slopes = compute(gldas, Time_periods);
	nldas_slopes = compute(nldas, Time_periods);
	grace_slopes = compute(grace, Time_periods);
	lsdm_slopes = compute(lsdm, Time_periods);
	gps_slopes = compute(gps, Time_periods);
	outputs(station, Time_periods, EQs, gldas, nldas, grace, lsdm, gps, gldas_slopes, nldas_slopes, grace_slopes, lsdm_slopes, gps_slopes);
	return;


def configure_times():
	T1="20060108:20100109"  # earthquake was 20050615. We want to skip six months
	T2="20100710:20140309"  # earthquake was 20100110. We want to skip six months. 
	T3="20140317:20161207"
	T4="20161215:20180915"
	Time_periods=[]; EQs=[]; 
	Time_periods.append(T1);
	Time_periods.append(T2);
	Time_periods.append(T3);
	Time_periods.append(T4);
	EQs.append(dt.datetime.strptime("20100110","%Y%m%d"));
	EQs.append(dt.datetime.strptime("20140310","%Y%m%d"));
	EQs.append(dt.datetime.strptime("20161208","%Y%m%d"));
	return Time_periods, EQs;

def inputs(station):
	[gldas,_,_] = gps_input_pipeline.get_station_data(station,'gldas');
	[nldas,_,_] = gps_input_pipeline.get_station_data(station,'nldas');
	[grace,_,_] = gps_input_pipeline.get_station_data(station,'grace');
	[lsdm,_,_] = gps_input_pipeline.get_station_data(station,'lsdm');
	[gps,offset_obj,EQtimes] = gps_input_pipeline.get_station_data(station,'pbo');
	gps=offsets.remove_offsets(gps, offset_obj);
	gps=offsets.remove_offsets(gps, EQtimes);
	gps=gps_seasonal_removals.make_detrended_ts(gps,0,'lssq');  # removing trend only
	return [gldas, nldas, grace, lsdm, gps];

def compute(ts_obj, Time_periods):
	slope_collection=[]; 
	if ts_obj==[]:
		return [];
	for i in range(len(Time_periods)):
		start_time=dt.datetime.strptime(Time_periods[i].split(":")[0],"%Y%m%d");
		end_time=dt.datetime.strptime(Time_periods[i].split(":")[1],"%Y%m%d");
		[east_slope, north_slope, vert_slope, east_std, north_std, vert_std] = gps_ts_functions.get_slope(ts_obj, starttime=start_time, endtime=end_time, missing_fraction=0.01);
		slope_collection.append(east_slope);
	return slope_collection; 


def outputs(station, Time_periods, EQs, gldas, nldas, grace, lsdm, gps, gldas_slopes, nldas_slopes, grace_slopes, lsdm_slopes, gps_slopes):
	start_time=dt.datetime.strptime("20060101","%Y%m%d");
	end_time = dt.datetime.strptime("20181201","%Y%m%d");
	offsets=[0,0,-1.5,1.5];  # one for each time period (for plotting only)
	markersize=0.9;
	f,axarr=plt.subplots(3,1,figsize=(14,10),sharex=True);

	# Dealing with error handling (you will end up with [] if there's been an error)
	plot_gldas=0; plot_nldas=0; plot_grace=0; plot_lsdm=0; 
	if gldas != []:
		plot_gldas=1; 
	if nldas != []:
		plot_nldas=1; 
	if grace != []:
		plot_grace=1;
	if lsdm != []:
		plot_lsdm=1; 

	# First plot: East component
	# Plotting that gets done no matter what. 
	axarr[0].plot_date(gps.dtarray, gps.dE, '.',color='gray',markersize=markersize,label='GPS');		
	axarr[0].set_xlim([start_time, end_time]);
	axarr[0].set_ylim([-5, 5]);
	axarr[0].set_ylabel('East (mm)');
	axarr[0].set_title(station+' loading models');
	for i in range(len(EQs)):
		lims = axarr[0].get_ylim();
		axarr[0].plot_date([EQs[i],EQs[i]],[lims[0], lims[1]],'g');

	if plot_gldas:
		axarr[0].plot_date(gldas.dtarray,gldas.dE,'--k',label='GLDAS');
		print("Plotting GLDAS");
	if plot_nldas:
		axarr[0].plot_date(nldas.dtarray,nldas.dE,'--r',label='NLDAS');
		print("Plotting NLDAS");
	if plot_grace:
		axarr[0].plot_date(grace.dtarray,grace.dE,'--b',label='GRACE');
		print("Plotting GRACE");
	if plot_lsdm:
		axarr[0].plot_date(lsdm.dtarray,lsdm.dE,'--g',label='LSDM');
		print("Plotting LSDM");
	

	# Slope annotations
	for i in range(len(Time_periods)):
		start=dt.datetime.strptime(Time_periods[i].split(":")[0],"%Y%m%d");
		end=dt.datetime.strptime(Time_periods[i].split(":")[1],"%Y%m%d");
		floats = gps_ts_functions.get_float_times([start, end]);		
		if plot_gldas:
			end_pos=offsets[i]+(floats[1]-floats[0])*gldas_slopes[i];
			axarr[0].plot([start, end],[offsets[i], end_pos],'k');		
		if plot_nldas:
			end_pos=offsets[i]+(floats[1]-floats[0])*nldas_slopes[i];
			axarr[0].plot([start, end],[offsets[i], end_pos],'r');
		if plot_grace:
			end_pos=offsets[i]+(floats[1]-floats[0])*grace_slopes[i];
			axarr[0].plot([start, end],[offsets[i], end_pos],'b');
		if plot_lsdm:
			end_pos=offsets[i]+(floats[1]-floats[0])*lsdm_slopes[i];
			axarr[0].plot([start, end],[offsets[i], end_pos],'g');		
		end_pos=offsets[i]+(floats[1]-floats[0])*gps_slopes[i];
		axarr[0].plot([start, end],[offsets[i], end_pos],color='gray');	
	
	for i in range(len(Time_periods)-1):
		end=dt.datetime.strptime(Time_periods[i].split(":")[1],"%Y%m%d");
		vertpos=2;
		axarr[0].text(end, vertpos+0.6, "%.3f mm/yr" % (gps_slopes[i+1]-gps_slopes[i]) ,color='gray');
		if plot_gldas: 
			axarr[0].text(end, vertpos, "%.3f mm/yr" % (gldas_slopes[i+1]-gldas_slopes[i]) ,color='k');
		if plot_nldas: 
			axarr[0].text(end, vertpos-0.6, "%.3f mm/yr" % (nldas_slopes[i+1]-nldas_slopes[i]) ,color='r');
		if plot_grace: 
			axarr[0].text(end, vertpos-1.2, "%.3f mm/yr" % (grace_slopes[i+1]-grace_slopes[i]) ,color='b');
		if plot_lsdm: 
			axarr[0].text(end, vertpos-1.8, "%.3f mm/yr" % (lsdm_slopes[i+1]-lsdm_slopes[i]) ,color='g');		
	axarr[0].legend();


	# Second plot
	axarr[1].plot_date(gps.dtarray, gps.dN, '.',color='gray',markersize=markersize,label='GPS');
	axarr[1].set_xlim([start_time, end_time]);
	axarr[1].set_ylim([-5, 5]);
	axarr[1].set_ylabel('North (mm)');
	for i in range(len(EQs)):
		lims = axarr[1].get_ylim();
		axarr[1].plot_date([EQs[i],EQs[i]],[lims[0], lims[1]],'g');	
	if plot_gldas:
		axarr[1].plot_date(gldas.dtarray,gldas.dN,'--k',label='GLDAS');
	if plot_nldas:
		axarr[1].plot_date(nldas.dtarray,nldas.dN,'--r',label='NLDAS');
	if plot_grace:
		axarr[1].plot_date(grace.dtarray,grace.dN,'--b',label='GRACE');
	if plot_lsdm:
		axarr[1].plot_date(lsdm.dtarray,lsdm.dN,'--g',label='LSDM');


	axarr[2].plot_date(gps.dtarray, gps.dU, '.',color='gray',markersize=markersize,label='GPS');
	axarr[2].set_xlim([start_time, end_time]);
	axarr[2].set_ylim([-20, 20]);
	axarr[2].set_ylabel('Up (mm)');
	for i in range(len(EQs)):
		lims = axarr[2].get_ylim();
		axarr[2].plot_date([EQs[i],EQs[i]],[lims[0], lims[1]],'g');	
	if plot_gldas:
		axarr[2].plot_date(gldas.dtarray,gldas.dU,'--k',label='GLDAS');
	if plot_nldas:
		axarr[2].plot_date(nldas.dtarray,nldas.dU,'--r',label='NLDAS');
	if plot_grace:
		axarr[2].plot_date(grace.dtarray,grace.dU,'--b',label='GRACE');
	if plot_lsdm:
		axarr[2].plot_date(lsdm.dtarray,lsdm.dU,'--g',label='LSDM');		

	plt.savefig("GRACEvsLDAS/"+station+"_model_compare.eps");
	return;


if __name__=="__main__":
	station="WDCB";
	compare_single_station(station);
