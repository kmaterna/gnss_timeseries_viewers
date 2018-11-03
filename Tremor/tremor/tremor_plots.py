# Here we make a cumulative plot with GPS
# Useful for viewing tremor. 

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import tremor_io
import tremor_tools


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
	[dt1, c1]=tremor_tools.get_cumulative_plot(tremor, box_interest1, start_time, end_time);
	[dt2, c2]=tremor_tools.get_cumulative_plot(tremor, box_interest2, start_time, end_time);
	[dt3, c3]=tremor_tools.get_cumulative_plot(tremor, box_interest3, start_time, end_time);

	station='P159';
	trend_out_gps=tremor_tools.get_detrended_gps_station(station);


	f,axarr=plt.subplots(2,1, sharex=True,figsize=(16,10));
	axarr[0].grid(True);
	axarr[0].plot_date(tremor.dtarray,tremor.latarray,'.',color='k',markersize=1);
	for item in eqtimes:
		axarr[0].plot_date([item, item],[tremor_latmin, tremor_latmax],color='red',linestyle='--',linewidth=2,marker=None);	
	axarr[0].set_xlim([start_time, end_time]);
	axarr[0].set_ylim([tremor_latmin, tremor_latmax]);
	axarr[0].set_ylabel('Latitude (degrees)',fontsize=20);
	axarr[0].tick_params(axis='both', which='major', labelsize=20);


	h1=axarr[1].plot_date(dt1,c1/max(c1),color='darkcyan',linestyle='-',linewidth=4,marker=None,label='coupling zone');
	h2=axarr[1].plot_date(dt2,c2/max(c2),color='darkorchid',linestyle='-',linewidth=4,marker=None,label='ETS zone');
	h3=axarr[1].plot_date(dt3,c3/max(c3),color='darkorange',linestyle='-',linewidth=4,marker=None,label='deep slip zone');
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


