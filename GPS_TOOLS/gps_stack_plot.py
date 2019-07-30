# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and eventually trend from the data. 
# Step 4: Plot in order of increasing latitude, colored by how close they are to the earthquake. 

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import collections
import sys
import datetime as dt 
import gps_io_functions
import gps_input_pipeline
import gps_ts_functions
import gps_seasonal_removals
import stations_within_radius
import offsets
import remove_ets_events


def driver():
	[stations, distances, EQtimes, proc_center, expname] = configure();
	[dataobj_list, offsetobj_list, eqobj_list] = inputs(stations, proc_center);
	[detrended_objects, stage1_objects, stage2_objects, sorted_distances, east_slope_obj] = compute(dataobj_list, offsetobj_list, eqobj_list, distances, EQtimes);
	
	output_full_ts(detrended_objects, sorted_distances, EQtimes, expname, "detrended", east_slope_obj);
	output_full_ts(stage1_objects, sorted_distances, EQtimes, expname, "noeq", east_slope_obj);
	output_full_ts(stage2_objects, sorted_distances, EQtimes, expname, "noeq_noseasons", east_slope_obj);
	vertical_plots(stage2_objects, sorted_distances, EQtimes, expname, "vertical");
	return;



def configure():
	EQcoords=[-125.134, 40.829]; expname='Mend'; radius = 120; 
	# EQcoords=[-123.0, 40.6]; expname='Mend'; radius = 75; 

	# EQcoords=[-124.0, 38.0];     expname='Nbay'; radius = 125; 
	# EQcoords=[-119.0, 34.5];     expname='SoCal';  radius = 25; # km
	# EQcoords=[-116.0, 34.5];     expname='Mojave';  radius = 35; # km
	# EQcoords=[-117.5, 35.5];     expname='ECSZ';  radius = 50; # km
	# EQcoords=[-117.5, 33.5];     expname='SD';  radius = 47; # km
	# EQcoords=[-119.0, 37.7];     expname='LVC';  radius = 30; # km
	EQtimes = [];  # What black lines do you want added to the figure? 
	EQtimes.append(dt.datetime.strptime("20140310", "%Y%m%d"));  # starts with the most important one
	EQtimes.append(dt.datetime.strptime("20050615", "%Y%m%d"));  # other earthquakes added to the figure
	EQtimes.append(dt.datetime.strptime("20100110", "%Y%m%d"));
	EQtimes.append(dt.datetime.strptime("20161208", "%Y%m%d"));
	
	proc_center='pbo';   # WHICH DATASTREAM DO YOU WANT?
 
	stations, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius);
	blacklist=["P316","P170","P158","TRND","P203","BBDM","KBRC","RYAN","BEAT"];
	stations_new=[]; distances_new=[];
	for i in range(len(stations)):
		if stations[i] not in blacklist:
			stations_new.append(stations[i]);
			distances_new.append(distances[i]);
	print(len(stations_new));
	print(stations_new);	
	return [stations_new, distances_new, EQtimes, proc_center, expname];


def inputs(station_names, proc_center):  # Returns a list of objects for time series data, offsets, and earthquakes
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[];
	for station_name in station_names:
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, proc_center, "NA");
		if myData.dtarray[-1]>dt.datetime.strptime("20140310","%Y%m%d") and myData.dtarray[0]<dt.datetime.strptime("20100310","%Y%m%d"):  
		# kicking out the stations that end early or start late. 
			dataobj_list.append(myData);
			offsetobj_list.append(offset_obj);
			eqobj_list.append(eq_obj);	
	return [dataobj_list, offsetobj_list, eqobj_list];


def compute(dataobj_list, offsetobj_list, eqobj_list, distances, EQtimes):

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_offsets = [x for _,x in sorted(zip(latitudes_list, offsetobj_list))];  # the raw, sorted data. 
	sorted_eqs = [x for _,x in sorted(zip(latitudes_list, eqobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

	# Detrended objects
	detrended_objects=[];
	for i in range(len(sorted_objects)):
		newobj=gps_seasonal_removals.make_detrended_ts(sorted_objects[i], 0, 'lssq');
		detrended_objects.append(newobj);

	# Objects with no earthquakes or seasonals
	stage1_objects = [];
	stage2_objects = [];
	east_slope_obj=[];
	for i in range(len(dataobj_list)):

		# Remove the steps earthquakes
		newobj=offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
		newobj=offsets.remove_offsets(newobj,sorted_eqs[i]);

		# The detrended TS without earthquakes
		stage1obj=gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq');
		stage1_objects.append(stage1obj);


		# NOTE: WRITTEN IN JUNE 2019
		# An experiment for removing ETS events
		# stage2obj=stage1obj;
		ets_intervals=remove_ets_events.input_tremor_days();
		stage2obj=gps_ts_functions.remove_outliers(stage1obj,3.0);  # 3 mm outlier def. 
		stage2obj=remove_ets_events.remove_ETS_times(stage2obj,ets_intervals, offset_num_days=15);  # 30 days on either end of the offsets
		stage2obj=gps_seasonal_removals.make_detrended_ts(stage2obj,0,'lssq');


		# The detrended TS without earthquakes or seasonals
		stage2obj=gps_seasonal_removals.make_detrended_ts(stage2obj, 1, 'lssq');

		stage2_objects.append(stage2obj);

		# Get the pre-event and post-event velocities (earthquakes removed)
		[east_slope_before, north_slope_before, vert_slope_before, _, _, _]=gps_ts_functions.get_slope(stage2obj,endtime=EQtimes[0]);
		[east_slope_after, north_slope_after, vert_slope_after, _, _, _]=gps_ts_functions.get_slope(stage2obj,starttime=EQtimes[0]);
		east_slope_after=np.round(east_slope_after,decimals=1);
		east_slope_before=np.round(east_slope_before,decimals=1);
		east_slope_obj.append([east_slope_before, east_slope_after]);

	return [detrended_objects, stage1_objects, stage2_objects, sorted_distances, east_slope_obj];




def output_full_ts(dataobj_list, distances, EQtimes, expname, filename, east_slope_obj):

	fig = plt.figure(figsize=(20,15),dpi=160);
	[f,axarr]=plt.subplots(1,2,sharex=True,sharey=True,figsize=(10,8))
	label_date=dt.datetime.strptime("20190215","%Y%m%d");
	start_time_plot=dt.datetime.strptime("20050101","%Y%m%d");
	end_time_plot=dt.datetime.strptime("20190116", "%Y%m%d");
	t1label=dt.datetime.strptime('20070701',"%Y%m%d");
	t2label=dt.datetime.strptime('20120101',"%Y%m%d");
	t3label=dt.datetime.strptime('20150501',"%Y%m%d");
	t4label=dt.datetime.strptime('20171201',"%Y%m%d");	

	offset=0;
	spacing=10;
	if expname=='Mend': # aesthetics only
		closest_station=70;
		farthest_station=120;
		# closest_station=10;
		# farthest_station=60;		
	else:
		closest_station=min(distances);  # km from event
		farthest_station=max(distances); # km from event
	color_boundary_object=matplotlib.colors.Normalize(vmin=closest_station,vmax=farthest_station, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# East
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		edata=dataobj_list[i].dE;
		edata=[x + offset for x in edata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = axarr[0].plot_date(dataobj_list[i].dtarray,edata,marker='+',markersize=1.5,color=line_color);
		axarr[0].text(label_date,offset,dataobj_list[i].name,fontsize=9,color=line_color);
		# axarr[0].text(dt.datetime.strptime("20050301", "%Y%m%d"),offset,east_slope_obj[i][0],fontsize=9,color='k');
		# axarr[0].text(EQtime,offset,east_slope_obj[i][1],fontsize=9,color='k');
	axarr[0].set_xlim(start_time_plot,end_time_plot);
	axarr[0].set_ylim([-10,offset+10])
	bottom,top=axarr[0].get_ylim();
	for i in range(len(EQtimes)):
		axarr[0].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	axarr[0].set_ylabel("East (mm)");
	axarr[0].set_title("East GPS Time Series")
	axarr[0].grid(True)
	axarr[0].text(t1label, top-spacing/2, "T1",fontsize=10, color='blue',fontweight='bold');
	axarr[0].text(t2label, top-spacing/2, "T2",fontsize=10, color='blue',fontweight='bold');
	axarr[0].text(t3label, top-spacing/2, "T3",fontsize=10, color='blue',fontweight='bold');
	axarr[0].text(t4label, top-spacing/2, "T4",fontsize=10, color='blue',fontweight='bold');	

	# North
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		ndata=dataobj_list[i].dN;
		ndata=[x + offset for x in ndata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = axarr[1].plot_date(dataobj_list[i].dtarray,ndata,marker='+',markersize=1.5, color=line_color);
	axarr[1].set_xlim(start_time_plot,end_time_plot);
	axarr[1].set_ylim([-10,offset+10])
	bottom,top=axarr[1].get_ylim();
	for i in range(len(EQtimes)):
		axarr[1].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	axarr[1].set_ylabel("North (mm)");
	axarr[1].set_title("North GPS Time Series")
	axarr[1].grid(True)
	axarr[1].text(t1label, top-spacing/2, "T1",fontsize=10, color='blue',fontweight='bold');
	axarr[1].text(t2label, top-spacing/2, "T2",fontsize=10, color='blue',fontweight='bold');
	axarr[1].text(t3label, top-spacing/2, "T3",fontsize=10, color='blue',fontweight='bold');
	axarr[1].text(t4label, top-spacing/2, "T4",fontsize=10, color='blue',fontweight='bold');	
	custom_cmap.set_array(range(int(closest_station),int(farthest_station)))
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilometers from 2014 Earthquake');

	plt.savefig(expname+'_Collective_TS_'+filename+'.jpg')	
	plt.close();
	print("Horizontal plots created.");

	return;


def vertical_plots(dataobj_list, distances, EQtimes, expname, filename):

	plt.figure(figsize=(6,8),dpi=160);
	label_date=dt.datetime.strptime("20190215","%Y%m%d");
	start_time_plot=dt.datetime.strptime("20050101","%Y%m%d");
	end_time_plot=dt.datetime.strptime("20190116", "%Y%m%d");

	offset=0;
	spacing=40;
	if expname=='Mend':  # aesthetics only
		closest_station=10;
		farthest_station=60;  # 70 to 120 was good for regular Mendocino
	else:
		closest_station=min(distances);  # km from event
		farthest_station=max(distances); # km from event
	color_boundary_object=matplotlib.colors.Normalize(vmin=closest_station,vmax=farthest_station, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# East
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		udata=dataobj_list[i].dU;
		udata=[x + offset for x in udata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = plt.gca().plot_date(dataobj_list[i].dtarray,udata,marker='+',markersize=1.5,color=line_color);
		plt.gca().text(label_date,offset,dataobj_list[i].name,fontsize=9,color=line_color);
		# axarr[0].text(dt.datetime.strptime("20050301", "%Y%m%d"),offset,east_slope_obj[i][0],fontsize=9,color='k');
		# axarr[0].text(EQtime,offset,east_slope_obj[i][1],fontsize=9,color='k');
	plt.gca().set_xlim(start_time_plot,end_time_plot);
	plt.gca().set_ylim([-10,offset+10])
	bottom,top=plt.gca().get_ylim();
	for i in range(len(EQtimes)):
		plt.gca().plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	plt.gca().set_ylabel("Vertical (mm)");
	plt.gca().set_title("Vertical GPS Time Series")
	plt.gca().grid(True)

	custom_cmap.set_array(range(int(closest_station),int(farthest_station)));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilometers from 2014 Earthquake');


	# new axis for plotting the map of california
	ax=plt.axes(position=[0.75,0.1,0.2,0.2],xticklabels=[],yticklabels=[]);
	[ca_lons,ca_lats]=np.loadtxt('../california_bdr',unpack=True);
	ax.plot(ca_lons,ca_lats,'k');
	for i in range(len(dataobj_list)):
		ax.plot(dataobj_list[i].coords[0],dataobj_list[i].coords[1],'.g',markersize=0.6);


	plt.savefig(expname+'_Collective_TS_'+filename+'.jpg');
	plt.close();
	print("Vertical plot created.");
	return;


