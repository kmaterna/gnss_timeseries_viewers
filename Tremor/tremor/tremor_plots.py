# Here we make a cumulative plot with GPS
# Useful for viewing tremor. 
# Added at the end: a second plot that has bins by depth. 

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import collections
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.dates as mdates
import tremor_io
import tremor_tools

# For reference:
TremorCat = collections.namedtuple("TremorCat",['dtarray','lonarray','latarray']);
TremorCatDepths = collections.namedtuple("TremorCatDepths",['dtarray','lonarray','latarray','depth']);

def simple_plot(tremor, tremortype):
	# Define bounds. 
	start_time=dt.datetime.strptime('20050301',"%Y%m%d");
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
	plt.savefig(tremortype+'_tremor_time_space.eps');
	return;


def timing_2014_plot(tremor, tremortype):
	# Looking at 2014 in more detail. 

	latmin=40;
	latmax=41;
	lonmin=-125;
	lonmax=-121.9;
	start_time=dt.datetime.strptime('20140308',"%Y%m%d");
	end_time=dt.datetime.strptime('20140316',"%Y%m%d");	
	eqtime14=dt.datetime.strptime('20140310 05',"%Y%m%d %H");  # time of 2014 earthquake
	eqtime16=dt.datetime.strptime('20161208 00',"%Y%m%d %H");  # time of 2016 earthquake

	# Make a restricted tremor catalog.
	box_dt=[]; box_lon=[]; box_lat=[];
	for i in range(len(tremor.dtarray)):
		if tremor.latarray[i]<latmax and tremor.latarray[i]>latmin:
			if tremor.dtarray[i]>=start_time and tremor.dtarray[i]<=end_time:
				box_dt.append(tremor.dtarray[i]);
				box_lon.append(tremor.lonarray[i]);
				box_lat.append(tremor.latarray[i]);

	# Cumulative plots
	box_interest1=[-124,-123.35,40,41];
	box_interest2=[-123.3,-122.8,40,41];
	box_interest3=[-122.8,-122,40,41];
	[dt1, c1]=tremor_tools.get_cumulative_plot(tremor, box_interest1, start_time, end_time);
	[dt2, c2]=tremor_tools.get_cumulative_plot(tremor, box_interest2, start_time, end_time);
	[dt3, c3]=tremor_tools.get_cumulative_plot(tremor, box_interest3, start_time, end_time);


	plt.figure(figsize=(10,10));
	plt.plot_date(box_dt,box_lon,'.',markersize=5);
	plt.plot_date([eqtime14,eqtime14],[lonmin,lonmax],'--k');
	plt.xlim([start_time, end_time]);
	plt.ylim([lonmin, lonmax]);
	plt.xlabel('Time');
	plt.ylabel('Longitude');
	plt.grid(True);
	
	# ax2=plt.gca().twinx();
	# h1=ax2.plot_date(dt1,c1/max(c1),color='darkcyan',linestyle='-',linewidth=4,marker=None,label='coupling zone');
	# h2=ax2.plot_date(dt2,c2/max(c2),color='darkorchid',linestyle='-',linewidth=4,marker=None,label='ETS zone');
	# h3=ax2.plot_date(dt3,c3/max(c3),color='darkorange',linestyle='-',linewidth=4,marker=None,label='deep slip zone');

	plt.savefig("2014_plot.eps");
	return;


def complex_plot(tremor,tremortype):
	# start_time=dt.datetime.strptime('20120301',"%Y%m%d");
	# end_time=dt.datetime.strptime('20181101',"%Y%m%d");
	start_time=dt.datetime.strptime('20060301',"%Y%m%d");
	end_time=dt.datetime.strptime('20141201',"%Y%m%d");	
	tremor_latmin=39;
	tremor_latmax=42.5;
	box_interest1=[-124,-123.35,40,41];
	box_interest2=[-123.3,-123,40,41];
	box_interest3=[-122.9,-122,40,41];
	eqtimes=[dt.datetime.strptime('20140310',"%Y%m%d"),
		dt.datetime.strptime('20161208',"%Y%m%d"),dt.datetime.strptime('20100110',"%Y%m%d")];

	# Cumulative plots. 
	[dt1, c1]=tremor_tools.get_cumulative_plot(tremor, box_interest1, start_time, end_time);
	[dt2, c2]=tremor_tools.get_cumulative_plot(tremor, box_interest2, start_time, end_time);
	[dt3, c3]=tremor_tools.get_cumulative_plot(tremor, box_interest3, start_time, end_time);

	station='P160';
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
	axarr[1].tick_params(axis='y', which='major', colors='black');
	axarr[1].tick_params(axis='both', which='major', labelsize=20);
	axarr[1].legend(loc=2,fontsize=18);
	plt.subplots_adjust(wspace=0, hspace=0.1)
	plt.savefig(tremortype+'_tremor_cumulative.eps');
	return;



def complex_plot_depths(tremor,tremortype):
	# This tremor object has depths associated. 
	# start_time=dt.datetime.strptime('20140201',"%Y%m%d");
	# end_time=dt.datetime.strptime('20140401',"%Y%m%d");  # the 2014 earthquake experiment
	# start_time=dt.datetime.strptime('20161101',"%Y%m%d");
	# end_time=dt.datetime.strptime('20170201',"%Y%m%d");  # the 2016 earthquake experiment
	# start_time=dt.datetime.strptime('20170725',"%Y%m%d");
	# end_time=dt.datetime.strptime('20170815',"%Y%m%d");  # the 2017 tremor experiment (shallow)
	# start_time=dt.datetime.strptime('20171215',"%Y%m%d");
	# end_time=dt.datetime.strptime('20180115',"%Y%m%d");  # the 2018 tremor experiment (shallow)
	start_time=dt.datetime.strptime('20120301',"%Y%m%d");  # the whole time period
	end_time=dt.datetime.strptime('20181101',"%Y%m%d");   # the whole time period
	# start_time=dt.datetime.strptime('20150101',"%Y%m%d");  # Aaron's new data
	# end_time=dt.datetime.strptime('20180101',"%Y%m%d");   # Aaron's new data

	# box_interest = [-125, -121, 40.1, 41];  # Nice
	box_interest = [-125, -121, 40.2, 40.8];  # Experiment
	# depth_interest1=[20, 24]; name1="20-24km";
	depth_interest2=[24, 35]; name2="24-35km";
	depth_interest3=[35, 65]; name3="35-65km";
	tremor_latmin=39;
	tremor_latmax=42.5;
	eqtimes=[dt.datetime.strptime('20140310',"%Y%m%d"),
		dt.datetime.strptime('20161208',"%Y%m%d"),dt.datetime.strptime('20100110',"%Y%m%d")];

	# Cumulative plots. 
	# [dt1, c1]=tremor_tools.get_cumulative_plot_depths(tremor, box_interest, depth_interest1, start_time, end_time);
	[dt2, c2]=tremor_tools.get_cumulative_plot_depths(tremor, box_interest, depth_interest2, start_time, end_time);
	[dt3, c3]=tremor_tools.get_cumulative_plot_depths(tremor, box_interest, depth_interest3, start_time, end_time);

	# Print the coordinates of the tremor in different depths, for GMT
	# shallowT=tremor_tools.restrict_to_box_depth(tremor, box_interest, depth_interest1, start_time, end_time);
	mediumT=tremor_tools.restrict_to_box_depth(tremor, box_interest, depth_interest2, start_time, end_time);
	deepT=tremor_tools.restrict_to_box_depth(tremor, box_interest, depth_interest3, start_time, end_time);
	# name1=name1+' (n='+str(len(shallowT.dtarray))+')';
	name2=name2+' (n='+str(len(mediumT.dtarray))+')';
	name3=name3+' (n='+str(len(deepT.dtarray))+')';
	# tremor_io.write_tremor_as_txt(shallowT, 'gmt/shallowrange.txt');
	tremor_io.write_tremor_as_txt(mediumT, 'gmt/medrange.txt');
	tremor_io.write_tremor_as_txt(deepT, 'gmt/deeprange.txt');

	station='P160';
	trend_out_gps=tremor_tools.get_detrended_gps_station(station);

	f,axarr=plt.subplots(2,1, sharex=True,figsize=(16,10));
	axarr[0].grid(True);
	axarr[0].plot_date(tremor.dtarray,tremor.latarray,'.',color='k',markersize=1);
	# axarr[0].plot_date(shallowT.dtarray,shallowT.latarray,'.',color='darkcyan',markersize=2);
	axarr[0].plot_date(mediumT.dtarray,mediumT.latarray,'.',color='darkorchid',markersize=1);
	for item in eqtimes:
		axarr[0].plot_date([item, item],[tremor_latmin, tremor_latmax],color='red',linestyle='--',linewidth=2,marker=None);	
	axarr[0].set_xlim([start_time, end_time]);
	axarr[0].set_ylim([tremor_latmin, tremor_latmax]);
	axarr[0].set_ylabel('Latitude (degrees)',fontsize=20);
	axarr[0].tick_params(axis='both', which='major', labelsize=20);


	# h1=axarr[1].plot_date(dt1,c1/max(c1),color='darkcyan',linestyle='-',linewidth=4,marker=None,label=name1);
	h2=axarr[1].plot_date(dt2,c2/max(c2),color='darkorchid',linestyle='-',linewidth=4,marker=None,label=name2);
	h3=axarr[1].plot_date(dt3,c3/max(c3),color='darkorange',linestyle='-',linewidth=4,marker=None,label=name3);
	axarr[1].text(dt.datetime.strptime("20130610","%Y%m%d"),1.02,'T2',color='red',fontsize=20,fontweight='bold');
	axarr[1].text(dt.datetime.strptime("20150610","%Y%m%d"),1.02,'T3',color='red',fontsize=20,fontweight='bold');
	axarr[1].text(dt.datetime.strptime("20171010","%Y%m%d"),1.02,'T4',color='red',fontsize=20,fontweight='bold');

	for item in eqtimes:
		axarr[1].plot_date([item, item],[0,max(c2)],color='red',linestyle='--',linewidth=2,marker=None);
	ax2=axarr[1].twinx();
	ax2.plot_date(trend_out_gps.dtarray, trend_out_gps.dE,marker='.',markersize=4,color='gray');
	ax2.tick_params(axis='both', which='major', labelsize=20);
	ax2.tick_params(axis='y', which='major', colors='gray');
	ax2.set_ylabel(station+' East (mm)',fontsize=20,color='gray');

	axarr[1].set_ylim([0,1]);
	axarr[1].set_ylabel('Norm. Tremor Counts',fontsize=20,color='black');
	axarr[1].grid(True);
	axarr[1].tick_params(axis='y', which='major', colors='black');
	axarr[1].tick_params(axis='both', which='major', labelsize=20);
	axarr[1].legend(loc=2,fontsize=18);
	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.savefig(tremortype+'_tremor_depth_cumulative.eps');
	return;


def histogram_depths(tremor, interval_list, box_interest, depth_interest):
	plt.figure();
	linecolor=[];
	# box_interest_hist=[-125, -121, 40.1, 40.9];  # a special histogram that focuses on the area near the coupling change patch
	# This wide one gives almost no difference between the ETS events. Beautiful plot. 
	box_interest_hist=[-125, -121, 40.2, 40.75];  # a special histogram that focuses on the area near the coupling change patch

	for i in range(len(interval_list)):
		start_time=interval_list[i][0];
		end_time=interval_list[i][1];
		tremor_box = tremor_tools.restrict_to_box_depth(tremor, box_interest_hist, depth_interest, start_time, end_time);
		a, b, c = plt.hist(tremor_box.depth,label=dt.datetime.strftime(start_time, "%Y-%m-%d"),histtype='step',density=True);
		linecolor.append(c[0].get_ec());
		plt.plot([np.mean(tremor_box.depth),np.median(tremor_box.depth)], [0, 1], color=linecolor[i], linestyle='--');
	
	plt.legend(loc='upper left');
	plt.ylim([0, 0.07]);
	plt.ylabel('Density');
	plt.xlabel('Depth (km)');
	plt.savefig('event_by_event/DepthHistogram.png');

	[ca_lon, ca_lat] = np.loadtxt("gmt/california_bdr",unpack=True);
	f, axarr = plt.subplots(2,4,figsize=(17,9));
	for i in range(len(interval_list)):
		if i<4:
			horiz_count=0; 
		else:
			horiz_count=1;
		start_time=interval_list[i][0];
		end_time=interval_list[i][1];
		tremor_box = tremor_tools.restrict_to_box_depth(tremor, box_interest, depth_interest, start_time, end_time);
		axarr[horiz_count][np.mod(i,4)].plot(ca_lon, ca_lat, color='k', linewidth=1);
		axarr[horiz_count][np.mod(i,4)].plot(tremor_box.lonarray, tremor_box.latarray,'s',label=dt.datetime.strftime(start_time, "%Y-%m-%d"),markersize=1.0,color=linecolor[i]);
		axarr[horiz_count][np.mod(i,4)].legend(loc='upper left');
		axarr[horiz_count][np.mod(i,4)].set_ylim([39.8, 41.2]);
		axarr[horiz_count][np.mod(i,4)].set_xlim([-124.4, -122.1]);
	plt.ylabel('Latitude');
	plt.xlabel('Longitude');
	plt.savefig('event_by_event/Map.png');
	return; 

def timespace_events(tremor, interval_list, box_interest, depth_interest):

	# Get the usual Python color palette
	color_count=[];
	testfit=plt.figure();
	for i in range(len(interval_list)):
		t1 = plt.plot([0,0],[0,0]);
		color_count.append(t1[0].get_color());

	f, axarr = plt.subplots(2,4,figsize=(17,9));
	eqtimes=[dt.datetime.strptime('20140310',"%Y%m%d"),
		dt.datetime.strptime('20161208',"%Y%m%d"),dt.datetime.strptime('20100110',"%Y%m%d")];
	[ca_lon, ca_lat] = np.loadtxt("gmt/california_bdr",unpack=True);
	ylim_min=39.8; 
	ylim_max=41.3;
	ylim_hist1=40.2;
	ylim_hist2=40.75;
	myFmt = mdates.DateFormatter('%b %d');

	for i in range(len(interval_list)):
		if i<4:
			horiz_count=0; 
		else:
			horiz_count=1;
		start_time=interval_list[i][0];
		end_time=interval_list[i][1];  # can hard-code the ETS length. 
		end_time_plot=start_time+dt.timedelta(days=75);
		tremor_box = tremor_tools.restrict_to_box_depth(tremor, box_interest, depth_interest, start_time, end_time_plot);
		axarr[horiz_count][np.mod(i,4)].fill_between([start_time, end_time],[ylim_min, ylim_min],[ylim_max,ylim_max],color='lightgray',alpha=0.7);
		axarr[horiz_count][np.mod(i,4)].plot_date(tremor_box.dtarray, tremor_box.latarray,'.',label=dt.datetime.strftime(start_time, "%Y-%m-%d"),markersize=1.0,color='k');
		axarr[horiz_count][np.mod(i,4)].legend(loc='upper left');
		for j in range(len(eqtimes)):
			axarr[horiz_count][np.mod(i,4)].plot_date([eqtimes[j],eqtimes[j]],[ylim_min, ylim_max],'--r');
		if np.mod(i,4)==0:
			axarr[horiz_count][np.mod(i,4)].set_ylabel('Longitude');
		axarr[horiz_count][np.mod(i,4)].plot_date([start_time, end_time],[ylim_hist1, ylim_hist1],'--k');
		axarr[horiz_count][np.mod(i,4)].plot_date([start_time, end_time],[ylim_hist2, ylim_hist2],'--k');
		axarr[horiz_count][np.mod(i,4)].xaxis.set_major_formatter(myFmt);
		plt.setp( axarr[horiz_count][np.mod(i,4)].xaxis.get_majorticklabels(), rotation=70, fontsize=9 );

		
		# California map with tremor
		inset_ax = inset_axes(axarr[horiz_count][np.mod(i,4)], height="33%", width="38%", loc=4);
		inset_ax.plot(tremor_box.lonarray, tremor_box.latarray,'s',label=dt.datetime.strftime(start_time, "%Y-%m-%d"),markersize=0.1,color=color_count[i]);
		inset_ax.plot(ca_lon, ca_lat, color='k', linewidth=1);
		inset_ax.set_ylim([39.8, 41.2]);
		inset_ax.set_xlim([-124.5, -122.1]);
		inset_ax.set_xticks([]);
		inset_ax.set_yticks([]);
		inset_ax.plot([-124.2, -122.1],[ylim_hist1, ylim_hist1],'--k',linewidth=0.5);
		inset_ax.plot([-124.2, -122.1],[ylim_hist2, ylim_hist2],'--k',linewidth=0.5);

		axarr[horiz_count][np.mod(i,4)].set_ylim([ylim_min, ylim_max]);
		axarr[horiz_count][np.mod(i,4)].set_xlim([start_time, end_time_plot]);
		axarr[horiz_count][np.mod(i,4)].grid(True);
	plt.savefig('event_by_event/timespace_events.png');
	plt.savefig('event_by_event/timespace_events.eps');
	return; 



if __name__=="__main__":
	tremortype='wech_custom';
	tremor=tremor_tools.read_custom_tremor(tremortype);
	complex_plot_depths(tremor, tremortype);
	# After this, you must go and make the GMT plots of the tremor (tremor_depth_ranges.sh)


