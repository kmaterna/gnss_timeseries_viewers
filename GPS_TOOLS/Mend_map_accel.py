# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers and earthquakes. Then identify 
# Step 4: Produce a table and plot of accelerations before/after time ranges. 

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
# Feature: verticals and horizontals at the same time, making two output plots
# Feature: feed seasonal type as parameter, and include that in the output_file name
# This lets us run several experiments.


import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import subprocess, sys
import gps_io_functions
import gps_ts_functions
import gps_seasonal_removals
import gps_input_pipeline
import offsets
import stations_within_radius
import haversine


def driver(EQcoords, outfile_name, deltat1, deltat2, fit_type, overall_size):
	[stations, map_coords, dt1_start, dt1_end, dt2_start, dt2_end, basename] = configure(EQcoords, outfile_name, deltat1, deltat2, fit_type, overall_size);
	[dataobj_list, offsetobj_list, eqobj_list] = inputs(stations);
	[noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj] = compute(dataobj_list, offsetobj_list, eqobj_list, dt1_start, dt1_end, dt2_start, dt2_end, fit_type);
	outputs(noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj, map_coords, basename);
	return;


def configure(EQcoords, outfile_name, deltat1, deltat2, fit_type, overall_size):
	basename=outfile_name+"_"+fit_type;
	dt1_start  = dt.datetime.strptime(deltat1[0], "%Y%m%d");
	dt1_end  = dt.datetime.strptime(deltat1[1], "%Y%m%d");
	dt2_start  = dt.datetime.strptime(deltat2[0], "%Y%m%d");
	dt2_end  = dt.datetime.strptime(deltat2[1], "%Y%m%d");

	if overall_size=='medium':
		radius=550;  # km.  
		map_coords=[EQcoords[0]-0.6, EQcoords[0]+6, EQcoords[1]-3.0, EQcoords[1]+3.0];
	elif overall_size=='huge':
		radius=-1;  # this is a special key for using a coordinate box instead of a radius
		map_coords=[-125.6, -110.0, 32.5, 48.5];
	else:
		map_coords=[EQcoords[0]-0.6, EQcoords[0]+4, EQcoords[1]-2.0, EQcoords[1]+2.0];
		radius=250;
	
	# Getting the stations of interest ('huge' means we just want within the box.)
	if radius==-1:
		stations = stations_within_radius.get_stations_within_box(map_coords);
	else:
		stations,_ = stations_within_radius.get_stations_within_radius(EQcoords, radius, map_coords);
	stations=gps_input_pipeline.remove_blacklist(stations);
	stations.append("CME6"); ## A special thing for CME6, not within PBO fields. 
	return [stations, map_coords, dt1_start, dt1_end, dt2_start, dt2_end, basename];

def inputs(station_names):
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[];
	for station_name in station_names:
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, 'pbo');
		dataobj_list.append(myData);
		offsetobj_list.append(offset_obj);
		eqobj_list.append(eq_obj);
	return [dataobj_list, offsetobj_list, eqobj_list];



def compute(dataobj_list, offsetobj_list, eqobj_list, dt1_start, dt1_end, dt2_start, dt2_end, fit_type):

	# No earthquakes objects
	noeq_objects = [];
	east_slope_obj=[];
	north_slope_obj=[];
	vert_slope_obj=[];
	period_after_start_date=7;  # wait a week. 

	# For the vertical correction. 
	names=[]; coords=[];
	for i in range(len(dataobj_list)):
		names.append(dataobj_list[i].name);
		coords.append(dataobj_list[i].coords);

	# The main processing loop for slopes. 
	for i in range(len(dataobj_list)):
		# Remove the earthquakes
		newobj=offsets.remove_antenna_offsets(dataobj_list[i],offsetobj_list[i]);
		newobj=offsets.remove_earthquakes(newobj, eqobj_list[i]);
		if fit_type=='none':
			newobj=gps_seasonal_removals.make_detrended_ts(newobj, 0, fit_type);  # remove seasonals
		else:
			newobj=gps_seasonal_removals.make_detrended_ts(newobj, 1, fit_type);  # remove seasonals

		# What is this? It looks like awful code. 
		if newobj.dN[0]==1.0000:
			print("Passing because we haven't computed GRACE yet...");
			noeq_objects.append(newobj);
			east_slope_obj.append([np.nan, np.nan]);
			north_slope_obj.append([np.nan, np.nan]);
			vert_slope_obj.append([np.nan, np.nan]);
			continue;

		noeq_objects.append(newobj);

		# Get the pre-event and post-event velocities (earthquakes removed)
		[east_slope_before, north_slope_before, vert_slope_before, esig0, nsig0, usig0]=gps_ts_functions.get_slope(newobj,starttime=dt1_start+dt.timedelta(days=period_after_start_date),endtime=dt1_end);
		[east_slope_after, north_slope_after, vert_slope_after, esig1, nsig1, usig1]=gps_ts_functions.get_slope(newobj,starttime=dt2_start+dt.timedelta(days=period_after_start_date),endtime=dt2_end);

		# When do we ignore stations? When their detrended time series have a large variance. 
		# print(dataobj_list[i].name);
		# print(esig0);
		# print(esig1);

		critical_value=5;  # mm/yr
		if abs(esig0)>critical_value or abs(nsig0)>critical_value or abs(esig1)>critical_value or abs(nsig1)>critical_value:
			print("Kicking station out...")
			print(dataobj_list[i].name);
			[east_slope_after, north_slope_after, vert_slope_after]=[np.nan,np.nan,np.nan];
			[east_slope_before, north_slope_before, vert_slope_before]=[np.nan,np.nan,np.nan];

		else:
			east_slope_after=np.round(east_slope_after,decimals=1);
			east_slope_before=np.round(east_slope_before,decimals=1);
			north_slope_after=np.round(north_slope_after,decimals=1);
			north_slope_before=np.round(north_slope_before,decimals=1);
			vert_slope_after=np.round(vert_slope_after,decimals=1);
			vert_slope_before=np.round(vert_slope_before,decimals=1);
		
		east_slope_obj.append([east_slope_before, east_slope_after]);
		north_slope_obj.append([north_slope_before, north_slope_after]);
		vert_slope_obj.append([vert_slope_before, vert_slope_after]); 

	# Adjusting verticals by a reference station. 
	vert_slope_obj = adjust_by_reference_stations(names, coords, vert_slope_obj);

	return [noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj];



def adjust_by_reference_stations(names, coords, slope_obj):
	# How do we adjust the verticals for large-scale drought signatures? 

	reference_station='P208';
	coord_box=[-123,-121,39,42];
	eq_coords=[-124.81, 40.53];
	radius=250;
	max_radius=350;
	reference_type='radius'  # options = 'radius','box','station'

	new_slope_obj=[];
	background_slopes_before=[];
	background_slopes_after =[];

	for i in range(len(names)):
		if reference_type=='station':
			if names[i]==reference_station:
				background_slopes_before.append(slope_obj[i][0]);
				background_slopes_after.append(slope_obj[i][1]);
		elif reference_type=='box':
			if coords[i][0]>coord_box[0] and coords[i][0]<coord_box[1]:
				if coords[i][1]>coord_box[2] and coords[i][1]<coord_box[3]:
					background_slopes_before.append(slope_obj[i][0]);
					background_slopes_after.append(slope_obj[i][1]);
		elif reference_type=='radius':
			mydistance=haversine.distance([coords[i][1],coords[i][0]],[eq_coords[1],eq_coords[0]]);
			if mydistance>radius and mydistance<max_radius:
				background_slopes_before.append(slope_obj[i][0]);
				background_slopes_after.append(slope_obj[i][1]);
	
	vert_reference_before=np.nanmean(background_slopes_before);
	vert_reference_after =np.nanmean(background_slopes_after);
	print("Vert slope before: %f " % vert_reference_before);
	print("Vert slope after: %f " % vert_reference_after);

	for i in range(len(slope_obj)):
		new_slope_obj.append([slope_obj[i][0]-vert_reference_before, slope_obj[i][1]-vert_reference_after]);
	
	return new_slope_obj;



def outputs(noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj, map_coords, basename):
	ofile1=open(basename+'.txt','w');
	for i in range(len(noeq_objects)):
		ofile1.write("%f %f %f %f 0 %f 0 0 %s\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1]-east_slope_obj[i][0], (north_slope_obj[i][1]-north_slope_obj[i][0]), vert_slope_obj[i][1]-vert_slope_obj[i][0], noeq_objects[i].name) );
	ofile1.close();
	subprocess.call(['./accel_map_gps.gmt',basename+'.txt',str(map_coords[0]),str(map_coords[1]),str(map_coords[2]),str(map_coords[3]),basename],shell=False);
	print('./accel_map_gps.gmt '+str(map_coords[0])+' '+str(map_coords[1])+' '+str(map_coords[2])+' '+str(map_coords[3])+' '+basename);
	return;

