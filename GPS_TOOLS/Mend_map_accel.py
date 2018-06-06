# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers and earthquakes. Then identify 
# Step 4: Produce a table and plot of accelerations before/after time0. 

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

# Config for the best parameters yet: 
# EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
# EQtime  = dt.datetime.strptime("20140310", "%Y%m%d");
# radius=280;  # km. 
# pre_event_duration = 2; # years
# post_event_duration = 2; # years
# map_coords=[EQcoords[0]-0.6, EQcoords[0]+3, EQcoords[1]-1.5, EQcoords[1]+1.5];

# Big Map of NorCal:
# EQcoords=[-123.834, 39.029]; # North Bay Area
# EQtime  = dt.datetime.strptime("20140310", "%Y%m%d");
# earthquakes_dir = earthquakes_dir="../GPS_POS_DATA/Event_Files/";
# offsets_dir = "../GPS_POS_DATA/Offsets/";
# radius=450;  # km. 
# pre_event_duration = 2; # years
# post_event_duration = 2; # years
# map_coords=[EQcoords[0]-0.6, EQcoords[0]+4, EQcoords[1]-2.5, EQcoords[1]+2.5];


import numpy as np 
import datetime as dt 
import subprocess
import gps_io_functions
import gps_ts_functions
import stations_within_radius


def driver():
	[stations, distances, filenames, EQtime, earthquakes_dir, offsets_dir, pre_event_duration, post_event_duration, map_coords,outfile_name] = configure();
	dataobj_list = inputs(stations, filenames);
	[noeq_objects, east_slope_obj, north_slope_obj] = compute(dataobj_list, distances, earthquakes_dir, offsets_dir, EQtime, pre_event_duration, post_event_duration);
	outputs(noeq_objects, east_slope_obj, north_slope_obj, map_coords,outfile_name);
	return;


def configure():
	# EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
	# EQcoords=[-122.834, 37.829]; # San Francisco Bay Area
	# EQcoords=[-123.834, 39.029]; # North Bay Area
	# EQcoords=[-125.134, 43.829]; # Oregon
	# EQcoords=[-125.134, 46.829]; # Washington
	EQcoords=[-125.134, 48.829]; # Canada
	EQtime  = dt.datetime.strptime("20140310", "%Y%m%d");
	earthquakes_dir = earthquakes_dir="../GPS_POS_DATA/Event_Files/";
	offsets_dir = "../GPS_POS_DATA/Offsets/";
	radius=450;  # km. 
	pre_event_duration = 2; # years
	post_event_duration = 2; # years
	map_coords=[EQcoords[0]-0.6, EQcoords[0]+4, EQcoords[1]-2.5, EQcoords[1]+2.5];
	outfile_name='accel_map.ps'
	stations, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius, map_coords);
	filenames=[];
	for station in stations:
		filenames.append("../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos");
	return [stations, distances, filenames, EQtime, earthquakes_dir, offsets_dir, pre_event_duration, post_event_duration, map_coords, outfile_name];




def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;

def compute(dataobj_list, distances, earthquakes_dir, offsets_dir, EQtime, pre_event_duration, post_event_duration):
	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

	# No earthquakes objects
	noeq_objects = [];
	east_slope_obj=[];
	north_slope_obj=[];
	for i in range(len(dataobj_list)):
		# Remove the earthquakes
		newobj=gps_ts_functions.remove_offsets(sorted_objects[i],offsets_dir);
		newobj=gps_ts_functions.remove_earthquakes(newobj,earthquakes_dir);
		noeq_objects.append(newobj);

		# Get the pre-event and post-event velocities (earthquakes removed)
		[east_slope_before, north_slope_before, vert_slope_before]=gps_ts_functions.get_slope(newobj,starttime=EQtime-dt.timedelta(days=pre_event_duration*365),endtime=EQtime);
		[east_slope_after, north_slope_after, vert_slope_after]=gps_ts_functions.get_slope(newobj,starttime=EQtime,endtime=EQtime+dt.timedelta(days=post_event_duration*365));
		east_slope_after=np.round(east_slope_after,decimals=1);
		east_slope_before=np.round(east_slope_before,decimals=1);
		east_slope_obj.append([east_slope_before, east_slope_after]);
		north_slope_after=np.round(north_slope_after,decimals=1);
		north_slope_before=np.round(north_slope_before,decimals=1);
		north_slope_obj.append([north_slope_before, north_slope_after]);
	return [noeq_objects, east_slope_obj, north_slope_obj];


def outputs(noeq_objects, east_slope_obj, north_slope_obj, map_coords,outfile_name):
	ofile=open('accelerations.txt','w');
	for i in range(len(noeq_objects)):
		ofile.write("%f %f %f %f 0 0 0\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1]-east_slope_obj[i][0], north_slope_obj[i][1]-north_slope_obj[i][0]) );
		# ofile.write("%f %f %f %f 0 0 0 %s\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1]-east_slope_obj[i][0], north_slope_obj[i][1]-north_slope_obj[i][0], noeq_objects[i].name) );
	ofile.close();
	subprocess.call(['./accel_map_gps.gmt',str(map_coords[0]),str(map_coords[1]),str(map_coords[2]),str(map_coords[3]),outfile_name],shell=False);
	print('./accel_map_gps.gmt '+str(map_coords[0])+' '+str(map_coords[1])+' '+str(map_coords[2])+' '+str(map_coords[3])+' '+outfile_name);
	return;

