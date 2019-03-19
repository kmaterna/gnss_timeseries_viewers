#   Goal: Compare 

import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import gps_input_pipeline
import gps_io_functions

def compare_single_station(station):
	[gldas, nldas, grace] = inputs(station);
	outputs(station, gldas, nldas, grace);
	return;

def inputs(station):
	[gldas,_,_] = gps_input_pipeline.get_station_data(station,'gldas');
	[nldas,_,_] = gps_input_pipeline.get_station_data(station,'nldas');
	[grace,_,_] = gps_input_pipeline.get_station_data(station,'grace');
	return [gldas, nldas, grace];

def compute_slopes():
	# Borrow this from the driver_east_plot.py
	return; 

def outputs(station, gldas, nldas, grace):
	start_time=dt.datetime.strptime("20120101","%Y%m%d");
	end_time = dt.datetime.strptime("20181201","%Y%m%d");
	f,axarr=plt.subplots(3,1,figsize=(14,8),sharex=True);

	axarr[0].plot_date(gldas.dtarray,gldas.dE,'--k',label='GLDAS');
	axarr[0].plot_date(nldas.dtarray,nldas.dE,'--r',label='NLDAS');
	axarr[0].plot_date(grace.dtarray,grace.dE,'--b',label='GRACE');
	axarr[0].set_xlim([start_time, end_time]);
	axarr[0].set_ylabel('East (mm)');
	axarr[0].legend();
	axarr[0].set_title(station+' loading models');
	
	axarr[1].plot_date(gldas.dtarray,gldas.dN,'--k',label='GLDAS');
	axarr[1].plot_date(nldas.dtarray,nldas.dN,'--r',label='NLDAS');
	axarr[1].plot_date(grace.dtarray,grace.dN,'--b',label='GRACE');
	axarr[1].set_xlim([start_time, end_time]);
	axarr[1].set_ylabel('North (mm)');
	axarr[1].legend();

	axarr[2].plot_date(gldas.dtarray,gldas.dU,'--k',label='GLDAS');
	axarr[2].plot_date(nldas.dtarray,nldas.dU,'--r',label='NLDAS');
	axarr[2].plot_date(grace.dtarray,grace.dU,'--b',label='GRACE');
	axarr[2].set_xlim([start_time, end_time]);
	axarr[2].set_ylabel('Up (mm)');
	axarr[2].legend();

	plt.savefig("GRACEvsLDAS/"+station+"_model_compare.eps");
	return;


if __name__=="__main__":
	station="P327";
	compare_single_station(station);
