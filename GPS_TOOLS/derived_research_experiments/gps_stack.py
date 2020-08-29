# Python viewing to see a stack of stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and eventually trend from the data. 
# Step 4: Plot in order of increasing latitude, colored by how close they are to the central point

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import matplotlib.pyplot as plt 
import collections
import subprocess
import sys
import datetime as dt 
import gps_input_pipeline
import gps_ts_functions
import gps_seasonal_removals
import stations_within_radius
import offsets
import remove_ets_events
import outputs_gps_stacks
import movie_tool


Parameters=collections.namedtuple("Parameters",['expname','proc_center','refframe','center','radius','stations','distances','blacklist','outdir', 'outname']);


def driver():
	myparams = configure();
	[dataobj_list, offsetobj_list, eqobj_list, paired_distances] = gps_input_pipeline.multi_station_inputs(myparams.stations, myparams.blacklist, myparams.proc_center, myparams.refframe, myparams.distances);
	[detrend_objects, no_offset_objects, no_offsets_no_trends, no_offsets_no_trends_no_seasons, sorted_distances] = compute(dataobj_list, offsetobj_list, eqobj_list, paired_distances);
	
	# outputs_gps_stacks.horizontal_full_ts(no_offsets_no_trends, sorted_distances, myparams, "noeq");
	outputs_gps_stacks.horizontal_full_ts(no_offsets_no_trends_no_seasons, sorted_distances, myparams, "noeq_noseasons");
	outputs_gps_stacks.vertical_full_ts(no_offsets_no_trends_no_seasons, sorted_distances, myparams);

	# outputs_gps_stacks.horizontal_filtered_plots(no_offsets_no_trends_no_seasons, sorted_distances, myparams);
	# outputs_gps_stacks.vertical_filtered_plots(no_offsets_no_trends_no_seasons, sorted_distances, myparams);
	# outputs_gps_stacks.vertical_filtered_plots(no_offset_objects, sorted_distances, myparams, 'trendsin_');
	# outputs_gps_stacks.pygmt_map(no_offsets_no_trends_no_seasons,myparams);

	movie_tool.movie_driver(no_offsets_no_trends_no_seasons, sorted_distances, myparams);  # make a movie (optional)

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
	center=[-115.5, 33]; expname='SSGF'; radius = 25; 
	# center=[-115.5, 33]; expname='SSGF'; radius =40; 

	proc_center='cwu';   # WHICH DATASTREAM DO YOU WANT?
	refframe = 'NA';     # WHICH REFERENCE FRAME? 

	stations, lons, lats, distances = stations_within_radius.get_stations_within_radius(center, radius, network=proc_center);
	blacklist=["P316","P170","P158","TRND","P203","BBDM","KBRC","RYAN","BEAT","CAEC","MEXI","BOMG","FSHB"];  # This is global, just keeps growing
	outdir=expname+"_"+proc_center
	subprocess.call(["mkdir","-p",outdir],shell=False);
	outname=expname+"_"+str(center[0])+"_"+str(center[1])+"_"+str(radius)
	myparams=Parameters(expname=expname, proc_center=proc_center, refframe=refframe, center=center, radius=radius, stations=stations, distances=distances, blacklist=blacklist, outdir=outdir, outname=outname);
	return myparams;


def compute(dataobj_list, offsetobj_list, eqobj_list, distances):

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_offsets = [x for _,x in sorted(zip(latitudes_list, offsetobj_list))];  # the raw, sorted data. 
	sorted_eqs = [x for _,x in sorted(zip(latitudes_list, eqobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

	detrended_objects = []; 
	no_offset_objects = [];
	no_offsets_no_trends = [];
	no_offsets_no_trends_no_seasons = [];

	# Detrended objects (or objects with trends and no offsets; depends on what you want.)
	for i in range(len(sorted_objects)):
		newobj=gps_seasonal_removals.make_detrended_ts(sorted_objects[i], 0, 'lssq');
		detrended_objects.append(newobj);  # still has offsets, doesn't have trends
		
		newobj=offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
		newobj=offsets.remove_offsets(newobj,sorted_eqs[i]);
		no_offset_objects.append(newobj);  # still has trends, doesn't have offsets

	# Objects with no earthquakes or seasonals
	for i in range(len(dataobj_list)):

		# Remove the steps earthquakes
		newobj=offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
		newobj=offsets.remove_offsets(newobj,sorted_eqs[i]);
		newobj=gps_ts_functions.remove_outliers(newobj, 20);  # 20mm outlier definition

		# The detrended TS without earthquakes
		stage1obj=gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq');
		no_offsets_no_trends.append(stage1obj);

		# The detrended TS without earthquakes or seasonals
		stage2obj=gps_seasonal_removals.make_detrended_ts(stage1obj, 1, 'lssq');
		no_offsets_no_trends_no_seasons.append(stage2obj);

	return [detrended_objects, no_offset_objects, no_offsets_no_trends, no_offsets_no_trends_no_seasons, sorted_distances];


		# # NOTE: WRITTEN IN JUNE 2019
		# # An experiment for removing ETS events
		# # stage2obj=stage1obj;
		# ets_intervals=remove_ets_events.input_tremor_days();
		# stage2obj=gps_ts_functions.remove_outliers(stage1obj,3.0);  # 3 mm outlier def. 
		# stage2obj=remove_ets_events.remove_ETS_times(stage2obj,ets_intervals, offset_num_days=15);  # 30 days on either end of the offsets
		# stage2obj=gps_seasonal_removals.make_detrended_ts(stage2obj,0,'lssq');


