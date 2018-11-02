# Here we make a cumulative plot with GPS
# Useful for viewing tremor. 

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import tremor_io
import gps_input_pipeline
import offsets
import gps_ts_functions
import gps_seasonal_removals


def simple_plot(tremor):
	# Define bounds. 
	start_time=dt.datetime.strptime('20120301',"%Y%m%d");
	end_time=dt.datetime.strptime('20181101',"%Y%m%d");
	latmin=39; latmax=42;

	# Make a very simple plot. 
	plt.figure(figsize=(16,10));
	plt.grid(True);
	plt.plot_date(tremor.dtarray,tremor.latarray,'.',color='k',markersize=1);
	plt.xlim([start_time, end_time]);
	plt.ylim([latmin, latmax]);
	plt.xlabel('Time',fontsize=20);
	plt.ylabel('Latitude (degrees)',fontsize=20);
	plt.tick_params(axis='both', which='major', labelsize=20);
	plt.savefig('tremor_time_space.eps');
	return;

def get_detrended_gps_station(station_name):
	datasource='pbo';
	[newData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, datasource);
	newData=offsets.remove_antenna_offsets(newData, offset_obj);
	newData=gps_ts_functions.remove_outliers(newData, 5);  # mm for outliers
	newData=offsets.remove_earthquakes(newData, eq_obj);
	trend_out=gps_seasonal_removals.make_detrended_ts(newData, 1, 'notch');
	return trend_out;


def get_cumulative_plot(tremor, box_interest, start_time, end_time):
	# Returns two arrays that can be plotted against each other to give the cumulative tremor plot. 
	dt_interest=[];
	cnumber=[];
	dt_interest.append(start_time);
	cnumber.append(0);
	for i in range(len(tremor.dtarray)):
		if tremor.dtarray[i]>start_time:
			if tremor.lonarray[i]>box_interest[0] and tremor.lonarray[i]<box_interest[1]:
				if tremor.latarray[i]>box_interest[2] and tremor.latarray[i]<box_interest[3]:
					dt_interest.append(tremor.dtarray[i]);
					dt_interest.append(tremor.dtarray[i]);
					cnumber.append(cnumber[-1]);
					cnumber.append(cnumber[-1]+1);
	cnumber=np.array(cnumber);	
	return [dt_interest, cnumber];


def complex_plot(tremor):
	start_time=dt.datetime.strptime('20120301',"%Y%m%d");
	end_time=dt.datetime.strptime('20181101',"%Y%m%d");
	tremor_latmin=39;
	tremor_latmax=42.5;
	box_interest1=[-124,-123.35,40,41];
	box_interest2=[-123.3,-123,40,41];
	box_interest3=[-122.9,-122,40,41];
	eqtimes=[dt.datetime.strptime('20140310',"%Y%m%d"),
		dt.datetime.strptime('20161208',"%Y%m%d")];

	# Cumulative plots. 
	[dt1, c1]=get_cumulative_plot(tremor, box_interest1, start_time, end_time);
	[dt2, c2]=get_cumulative_plot(tremor, box_interest2, start_time, end_time);
	[dt3, c3]=get_cumulative_plot(tremor, box_interest3, start_time, end_time);

	station='P159';
	trend_out_gps=get_detrended_gps_station(station);


	f,axarr=plt.subplots(2,1, sharex=True,figsize=(16,10));
	axarr[0].grid(True);
	axarr[0].plot_date(tremor.dtarray,tremor.latarray,'.',color='k',markersize=1);
	for item in eqtimes:
		axarr[0].plot_date([item, item],[tremor_latmin, tremor_latmax],color='red',linestyle='--',linewidth=2,marker=None);	
	axarr[0].set_xlim([start_time, end_time]);
	axarr[0].set_ylim([tremor_latmin, tremor_latmax]);
	axarr[0].set_ylabel('Latitude (degrees)',fontsize=20);
	axarr[0].tick_params(axis='both', which='major', labelsize=20);


	h1=axarr[1].plot_date(dt1,c1/max(c1),color='darkcyan',linestyle='-',linewidth=4,marker=None,label='18-27km (flat)');
	h2=axarr[1].plot_date(dt2,c2/max(c2),color='darkorchid',linestyle='-',linewidth=4,marker=None,label='28-38km (steep)');
	h3=axarr[1].plot_date(dt3,c3/max(c3),color='darkorange',linestyle='-',linewidth=4,marker=None,label='>40km (steep)');
	for item in eqtimes:
		axarr[1].plot_date([item, item],[0,max(c1)],color='red',linestyle='--',linewidth=2,marker=None);
	ax2=axarr[1].twinx();
	ax2.plot_date(trend_out_gps.dtarray, trend_out_gps.dE,marker='.',markersize=4,color='gray');
	ax2.tick_params(axis='both', which='major', labelsize=20);
	ax2.tick_params(axis='y', which='major', colors='gray');
	ax2.set_ylabel(station+' East (mm)',fontsize=20,color='gray');

	axarr[1].set_ylim([0,1]);
	axarr[1].set_ylabel('Norm. Tremor Counts',fontsize=20,color='black');
	axarr[1].grid(True);
	axarr[1].set_xlabel('Time',fontsize=20);
	axarr[1].tick_params(axis='y', which='major', colors='black');
	axarr[1].tick_params(axis='both', which='major', labelsize=20);
	axarr[1].legend(loc=2,fontsize=18);
	plt.subplots_adjust(wspace=0, hspace=0.1)
	plt.savefig('tremor_cumulative.eps');
	return;



if __name__=="__main__":
	tremor=tremor_io.read_wech("../../GPS_POS_DATA/tremor/08_01_2009_10_31_2018.txt");
	complex_plot(tremor);


