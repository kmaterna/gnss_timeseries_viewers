# Python viewing to see a stack of stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and eventually trend from the data. 
# Step 4: Plot in order of increasing latitude, colored by how close they are to the central point

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import scipy.ndimage
import collections
import subprocess
import sys
import datetime as dt 
import gps_io_functions
import gps_input_pipeline
import gps_ts_functions
import gps_seasonal_removals
import stations_within_radius
import offsets
import remove_ets_events
import pygmt


Parameters=collections.namedtuple("Parameters",['expname','proc_center','center','radius','stations','distances','blacklist','outdir', 'outname']);


def driver():
	myparams = configure();
	[dataobj_list, offsetobj_list, eqobj_list, paired_distances] = gps_input_pipeline.multi_station_inputs(myparams.stations, myparams.blacklist, myparams.proc_center, myparams.distances);
	[detrend_objects, stage1_objects, stage2_objects, sorted_distances] = compute(dataobj_list, offsetobj_list, eqobj_list, paired_distances);
	output_full_ts(detrend_objects, sorted_distances, myparams, "detrended");
	output_full_ts(stage1_objects, sorted_distances, myparams, "noeq");
	output_full_ts(stage2_objects, sorted_distances, myparams, "noeq_noseasons");
	vertical_plots(stage2_objects, sorted_distances, myparams);
	vertical_filtered_plots(stage2_objects, sorted_distances, myparams);
	# vertical_filtered_plots(detrend_objects, sorted_distances, myparams);
	pygmt_map(stage2_objects,myparams);
	return;



def configure():
	# center=[-125.134, 40.829]; expname='Mend'; radius = 120; # Keeper for the 2014 earthquake
	# center=[-125.134, 40.829]; expname='Mend'; radius = 90; # 
	# center=[-124, 40.5]; expname='Humboldt'; radius = 60; # 
	# center=[-122.5, 40.5]; expname='Chico'; radius = 75; # 
	# center=[-124.0, 38.0];     expname='Nbay'; radius = 125; 
	# center=[-119.0, 34.5];     expname='SoCal';  radius = 25; # km
	# center=[-116.0, 34.5];     expname='Mojave';  radius = 35; # km
	# center=[-117.5, 35.5];     expname='ECSZ';  radius = 50; # km
	# center=[-119.0, 37.7];     expname='LVC';  radius = 30; # km
	# center=[-115.5, 32.85]; expname='SSGF'; radius = 20; 
	center=[-115.5, 33]; expname='SSGF'; radius = 80; 

	proc_center='cwu';   # WHICH DATASTREAM DO YOU WANT?

	stations, distances = stations_within_radius.get_stations_within_radius(center, radius, network=proc_center);
	blacklist=["P316","P170","P158","TRND","P203","BBDM","KBRC","RYAN","BEAT","CAEC","MEXI","BOMG"];  # This is global, just keeps growing
	outdir=expname+"_"+proc_center
	subprocess.call(["mkdir","-p",outdir],shell=False);
	outname=expname+"_"+str(center[0])+"_"+str(center[1])+"_"+str(radius)
	myparams=Parameters(expname=expname, proc_center=proc_center, center=center, radius=radius, stations=stations, distances=distances, blacklist=blacklist, outdir=outdir, outname=outname);
	return myparams;


def compute(dataobj_list, offsetobj_list, eqobj_list, distances):

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_offsets = [x for _,x in sorted(zip(latitudes_list, offsetobj_list))];  # the raw, sorted data. 
	sorted_eqs = [x for _,x in sorted(zip(latitudes_list, eqobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

	# Detrended objects
	detrended_objects=[];
	for i in range(len(sorted_objects)):
		newobj=gps_seasonal_removals.make_detrended_ts(sorted_objects[i], 0, 'lssq');
		# detrended_objects.append(newobj);
		detrended_objects.append(sorted_objects[i]);  # if you want trends in your vertical data

	# Objects with no earthquakes or seasonals
	stage1_objects = [];
	stage2_objects = [];
	for i in range(len(dataobj_list)):

		# Remove the steps earthquakes
		newobj=offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
		newobj=offsets.remove_offsets(newobj,sorted_eqs[i]);

		# The detrended TS without earthquakes
		stage1obj=gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq');
		stage1_objects.append(stage1obj);

		# The detrended TS without earthquakes or seasonals
		stage2obj=gps_seasonal_removals.make_detrended_ts(stage1obj, 1, 'lssq');
		stage2_objects.append(stage2obj);

	return [detrended_objects, stage1_objects, stage2_objects, sorted_distances];


		# # NOTE: WRITTEN IN JUNE 2019
		# # An experiment for removing ETS events
		# # stage2obj=stage1obj;
		# ets_intervals=remove_ets_events.input_tremor_days();
		# stage2obj=gps_ts_functions.remove_outliers(stage1obj,3.0);  # 3 mm outlier def. 
		# stage2obj=remove_ets_events.remove_ETS_times(stage2obj,ets_intervals, offset_num_days=15);  # 30 days on either end of the offsets
		# stage2obj=gps_seasonal_removals.make_detrended_ts(stage2obj,0,'lssq');




def configure_beautiful_plots(expname, distances):
	
	EQtimes=[]; labeltimes=[]; labels=[];

	if expname=='Mend' or expname=="Humboldt":
	# What black lines do you want added to the figure? This is good for Mendocino
		EQtimes.append(dt.datetime.strptime("20140310", "%Y%m%d"));  # starts with the most important one
		EQtimes.append(dt.datetime.strptime("20050615", "%Y%m%d"));  # other earthquakes added to the figure
		EQtimes.append(dt.datetime.strptime("20100110", "%Y%m%d"));
		EQtimes.append(dt.datetime.strptime("20161208", "%Y%m%d"));

	elif expname=="SSGF":
		# # This is good for SoCal
		EQtimes.append(dt.datetime.strptime("20100403", "%Y%m%d"));  # starts with the most important one
		EQtimes.append(dt.datetime.strptime("20050615", "%Y%m%d"));  # other earthquakes added to the figure
		EQtimes.append(dt.datetime.strptime("20120808", "%Y%m%d"));	

	if expname=="Mend":
		closest_station=70;
		farthest_station=120;
	else:
		closest_station=min(distances);  # km from event
		farthest_station=max(distances); # km from event	

	# Labels
	if expname=="Mend" or expname=="Humboldt":
		labeltimes.append(dt.datetime.strptime('20070701',"%Y%m%d"));
		labeltimes.append(dt.datetime.strptime('20120101',"%Y%m%d"));
		labeltimes.append(dt.datetime.strptime('20150501',"%Y%m%d"));
		labeltimes.append(dt.datetime.strptime('20171201',"%Y%m%d"));
		labels.append("T1");
		labels.append("T2");
		labels.append("T3");
		labels.append("T4");

	return EQtimes, labeltimes, labels, closest_station, farthest_station;



def output_full_ts(dataobj_list, distances, myparams, filename):

	fig = plt.figure(figsize=(20,15),dpi=160);
	[f,axarr]=plt.subplots(1,2,sharex=True,sharey=True,figsize=(10,8))
	label_date=dt.datetime.strptime("20200215","%Y%m%d");
	start_time_plot=dt.datetime.strptime("20050101","%Y%m%d");
	end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d");


	offset=0;
	spacing=10;
	EQtimes, labeltimes, labels, closest_station, farthest_station=configure_beautiful_plots(myparams.expname, distances);

	color_boundary_object=matplotlib.colors.Normalize(vmin=closest_station,vmax=farthest_station, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# East
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		edata=dataobj_list[i].dE;
		emean=np.nanmean(dataobj_list[i].dE);
		edata=[x + offset - emean for x in edata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = axarr[0].plot_date(dataobj_list[i].dtarray,edata,marker='+',markersize=1.5,color=line_color);
		axarr[0].text(label_date,offset,dataobj_list[i].name,fontsize=9,color=line_color);
	axarr[0].set_xlim(start_time_plot,end_time_plot);
	axarr[0].set_ylim([-30,offset+20])
	bottom,top=axarr[0].get_ylim();
	for i in range(len(EQtimes)):
		axarr[0].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	for i in range(len(labeltimes)):
		axarr[0].text(labeltimes[i], top-spacing/2, labels[i],fontsize=10, color='blue',fontweight='bold');	
	axarr[0].set_ylabel("East (mm)");
	axarr[0].set_title("East GPS Time Series")
	axarr[0].grid(True)

	# North
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		ndata=dataobj_list[i].dN;
		nmean=np.nanmean(dataobj_list[i].dN);
		ndata=[x + offset - nmean for x in ndata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = axarr[1].plot_date(dataobj_list[i].dtarray,ndata,marker='+',markersize=1.5, color=line_color);
	axarr[1].set_xlim(start_time_plot,end_time_plot);
	axarr[1].set_ylim([-30,offset+20])
	bottom,top=axarr[1].get_ylim();
	for i in range(len(EQtimes)):
		axarr[1].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	axarr[1].set_ylabel("North (mm)");
	axarr[1].set_title("North GPS Time Series")
	axarr[1].grid(True)
	for i in range(len(labeltimes)):
		axarr[1].text(labeltimes[i], top-spacing/2, labels[i],fontsize=10, color='blue',fontweight='bold');	
	custom_cmap.set_array(range(int(closest_station),int(farthest_station)))
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilometers from center');

	plt.savefig(myparams.outdir+"/"+myparams.outname+'_TS_'+filename+'.jpg')
	plt.close();
	print("Horizontal plots created.");

	return;


def vertical_plots(dataobj_list, distances, myparams):

	plt.figure(figsize=(6,8),dpi=160);
	label_date=dt.datetime.strptime("20200215","%Y%m%d");
	start_time_plot=dt.datetime.strptime("20050101","%Y%m%d");
	end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d");

	offset=0;
	spacing=40;
	EQtimes, labeltimes, labels, closest_station, farthest_station=configure_beautiful_plots(myparams.expname, distances);
	color_boundary_object=matplotlib.colors.Normalize(vmin=closest_station,vmax=farthest_station, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# Vertical
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		udata=dataobj_list[i].dU;
		umean=np.nanmean(dataobj_list[i].dU)
		udata=[x + offset - umean for x in udata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = plt.gca().plot_date(dataobj_list[i].dtarray,udata,marker='+',markersize=1.5,color=line_color);
		plt.gca().text(label_date,offset,dataobj_list[i].name,fontsize=9,color=line_color);
	plt.gca().set_xlim(start_time_plot,end_time_plot);
	plt.gca().set_ylim([-20,offset+20])
	bottom,top=plt.gca().get_ylim();
	for i in range(len(EQtimes)):
		plt.gca().plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	plt.gca().set_ylabel("Vertical (mm)");
	plt.gca().set_title("Vertical GPS Time Series")
	plt.gca().grid(True)

	custom_cmap.set_array(range(int(closest_station),int(farthest_station)));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilometers from center');

	# new axis for plotting the map of california
	ax=plt.axes(position=[0.8,0.1,0.2,0.2],xticklabels=[],yticklabels=[]);
	[ca_lons,ca_lats]=np.loadtxt('../california_bdr',unpack=True);
	ax.plot(ca_lons,ca_lats,'k');
	for i in range(len(dataobj_list)):
		ax.plot(dataobj_list[i].coords[0],dataobj_list[i].coords[1],'.g',markersize=0.6);

	# new axis for extra labels
	ax=plt.axes(position=[0.8,0.4,0.2,0.2],xticklabels=[],yticklabels=[]);
	ax.set_ylim([-1,1])
	ax.set_xlim([-0.1,1])
	ax.text(0,0.75,myparams.proc_center+" data");
	ax.text(0,0.37,myparams.center);
	ax.text(0,0,str(myparams.radius)+" km radius");

	plt.savefig(myparams.outdir+"/"+myparams.outname+'_TS_'+"vertical"+'.jpg');
	plt.close();
	print("Vertical plot created.");
	return;


def vertical_filtered_plots(dataobj_list, distances, myparams):

	plt.figure(figsize=(15,8),dpi=160);
	label_date=dt.datetime.strptime("20200215","%Y%m%d");
	start_time_plot=dt.datetime.strptime("20050101","%Y%m%d");
	end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d");

	spacing=40;
	EQtimes, labeltimes, labels, closest_station, farthest_station=configure_beautiful_plots(myparams.expname, distances);
	color_boundary_object=matplotlib.colors.Normalize(vmin=closest_station,vmax=farthest_station, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# Vertical
	for i in range(len(dataobj_list)):
		umean=np.mean(dataobj_list[i].dU);  # start at the mean. 
		# umean=dataobj_list[i].dU[0];  # start at the beginning
		line_color=custom_cmap.to_rgba(distances[i]);
		# l1 = plt.gca().plot(dataobj_list[i].dtarray,dataobj_list[i].dU-umean,linestyle='solid',linewidth=0,marker='.',color='red' ); # for debugging the filter
		udata=scipy.ndimage.median_filter(dataobj_list[i].dU-umean,size=365);
		l1 = plt.gca().plot(dataobj_list[i].dtarray,udata,linestyle='solid',linewidth=1,color=line_color );
		# plt.gca().text(label_date,offset,dataobj_list[i].name,fontsize=9,color=line_color);
	plt.gca().set_xlim(start_time_plot,end_time_plot);
	bottom,top=plt.gca().get_ylim();
	for i in range(len(EQtimes)):
		plt.gca().plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	plt.gca().set_ylabel("Filtered Vertical (mm)");
	plt.gca().set_title("Filtered Vertical GPS Time Series")
	plt.gca().grid(True)

	custom_cmap.set_array(range(int(closest_station),int(farthest_station)));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilometers from center');

	# new axis for plotting the map of california
	ax=plt.axes(position=[0.8,0.1,0.1,0.2],xticklabels=[],yticklabels=[]);
	[ca_lons,ca_lats]=np.loadtxt('../california_bdr',unpack=True);
	ax.plot(ca_lons,ca_lats,'k');
	for i in range(len(dataobj_list)):
		ax.plot(dataobj_list[i].coords[0],dataobj_list[i].coords[1],'.g',markersize=0.6);

	# new axis for extra labels
	ax=plt.axes(position=[0.8,0.4,0.1,0.2],xticklabels=[],yticklabels=[]);
	ax.set_ylim([-1,1])
	ax.set_xlim([-0.1,1])
	ax.text(0,0.75,myparams.proc_center+" data");
	ax.text(0,0.37,myparams.center);
	ax.text(0,0,str(myparams.radius)+" km radius");

	plt.savefig(myparams.outdir+"/"+myparams.outname+'_TS_'+"vertical_filt"+'.jpg');
	plt.close();
	print("Vertical plot created.");
	return;


def pygmt_map(ts_objects,myparams):

	offset=0.2;

	lons=[]; lats=[]; names=[];
	for i in range(len(ts_objects)):
		lons.append(ts_objects[i].coords[0]);
		lats.append(ts_objects[i].coords[1]);
		names.append(ts_objects[i].name);
	region=[min(lons)-offset,max(lons)+offset,min(lats)-offset,max(lats)+offset];		

	fig = pygmt.Figure()
	fig.basemap(region=region,projection="M8i",B="0.25");
	# fig.grdimage("@earth_relief_30s",region=region,I="+d");  # takes a little while the first time, but faster each time afterwards
	fig.coast(shorelines="0.5p,black",G='peachpuff2',S='skyblue',D="h");
	fig.coast(N='1',W='1.0p,black');
	fig.coast(N='2',W='0.5p,black');
	fig.text(x=[i+0.035 for i in lons],y=lats,text=names,font='15p,Helvetica-Bold,black');
	fig.plot(x=lons,y=lats,S='c0.1i',G='black',W='0.5p,black')
	fig.plot(x=myparams.center[0],y=myparams.center[1],S='a0.1i',G='red',W='0.5p,red')
	fig.savefig(myparams.outdir+"/"+myparams.outname+'_map.png');
	return;


