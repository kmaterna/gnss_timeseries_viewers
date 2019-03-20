#   Goal: Compare 

import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import gps_input_pipeline
import gps_io_functions
import gps_ts_functions

def compare_single_station(station):
	Time_periods, EQs=configure_times();
	[gldas, nldas, grace] = inputs(station);
	gldas_slopes = compute(gldas, Time_periods);
	nldas_slopes = compute(nldas, Time_periods);
	grace_slopes = compute(grace, Time_periods);
	outputs(station, Time_periods, EQs, gldas, nldas, grace, gldas_slopes, nldas_slopes, grace_slopes);
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
	return [gldas, nldas, grace];


def compute(ts_obj, Time_periods):
	slope_collection=[]; 
	for i in range(len(Time_periods)):
		start_time=dt.datetime.strptime(Time_periods[i].split(":")[0],"%Y%m%d");
		end_time=dt.datetime.strptime(Time_periods[i].split(":")[1],"%Y%m%d");
		[east_slope, north_slope, vert_slope, east_std, north_std, vert_std] = gps_ts_functions.get_slope(ts_obj, starttime=start_time, endtime=end_time, missing_fraction=0.01);
		slope_collection.append(east_slope);
	return slope_collection; 


def outputs(station, Time_periods, EQs, gldas, nldas, grace, gldas_slopes, nldas_slopes, grace_slopes):
	start_time=dt.datetime.strptime("20120101","%Y%m%d");
	end_time = dt.datetime.strptime("20181201","%Y%m%d");
	offsets=[0,0,-1.5,1.5];  # one for each time period (for plotting only)
	f,axarr=plt.subplots(3,1,figsize=(14,10),sharex=True);

	axarr[0].plot_date(gldas.dtarray,gldas.dE,'--k',label='GLDAS');
	axarr[0].plot_date(nldas.dtarray,nldas.dE,'--r',label='NLDAS');
	axarr[0].plot_date(grace.dtarray,grace.dE,'--b',label='GRACE');
	axarr[0].set_xlim([start_time, end_time]);
	axarr[0].set_ylabel('East (mm)');
	# Slope annotations
	for i in range(len(Time_periods)):
		start=dt.datetime.strptime(Time_periods[i].split(":")[0],"%Y%m%d");
		end=dt.datetime.strptime(Time_periods[i].split(":")[1],"%Y%m%d");
		floats = gps_ts_functions.get_float_times([start, end]);
		end_pos=offsets[i]+(floats[1]-floats[0])*gldas_slopes[i];
		axarr[0].plot([start, end],[offsets[i], end_pos],'k');
		end_pos=offsets[i]+(floats[1]-floats[0])*nldas_slopes[i];
		axarr[0].plot([start, end],[offsets[i], end_pos],'r');
		end_pos=offsets[i]+(floats[1]-floats[0])*grace_slopes[i];
		axarr[0].plot([start, end],[offsets[i], end_pos],'b');
	for i in range(len(EQs)):
		lims = axarr[0].get_ylim();
		axarr[0].plot_date([EQs[i],EQs[i]],[lims[0], lims[1]],'g');
	for i in range(len(Time_periods)-1):
		end=dt.datetime.strptime(Time_periods[i].split(":")[1],"%Y%m%d");
		vertpos=2;
		axarr[0].text(end, vertpos, "%.3f mm/yr" % (gldas_slopes[i+1]-gldas_slopes[i]) ,color='k');
		axarr[0].text(end, vertpos-0.3, "%.3f mm/yr" % (nldas_slopes[i+1]-nldas_slopes[i]) ,color='r');
		axarr[0].text(end, vertpos-0.6, "%.3f mm/yr" % (grace_slopes[i+1]-grace_slopes[i]) ,color='b');
	axarr[0].legend();
	axarr[0].set_title(station+' loading models');

	axarr[1].plot_date(gldas.dtarray,gldas.dN,'--k',label='GLDAS');
	axarr[1].plot_date(nldas.dtarray,nldas.dN,'--r',label='NLDAS');
	axarr[1].plot_date(grace.dtarray,grace.dN,'--b',label='GRACE');
	axarr[1].set_xlim([start_time, end_time]);
	axarr[1].set_ylabel('North (mm)');
	for i in range(len(EQs)):
		lims = axarr[1].get_ylim();
		axarr[1].plot_date([EQs[i],EQs[i]],[lims[0], lims[1]],'g');	


	axarr[2].plot_date(gldas.dtarray,gldas.dU,'--k',label='GLDAS');
	axarr[2].plot_date(nldas.dtarray,nldas.dU,'--r',label='NLDAS');
	axarr[2].plot_date(grace.dtarray,grace.dU,'--b',label='GRACE');
	axarr[2].set_xlim([start_time, end_time]);
	axarr[2].set_ylabel('Up (mm)');
	for i in range(len(EQs)):
		lims = axarr[2].get_ylim();
		axarr[2].plot_date([EQs[i],EQs[i]],[lims[0], lims[1]],'g');	

	plt.savefig("GRACEvsLDAS/"+station+"_model_compare.eps");
	return;


if __name__=="__main__":
	station="P332";
	compare_single_station(station);
