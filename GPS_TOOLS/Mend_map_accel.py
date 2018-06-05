# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers and earthquakes. Then identify 
# Step 4: Produce a table and plot of accelerations before/after time0. 

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import datetime as dt 
import subprocess
import gps_io_functions
import gps_ts_functions
import stations_within_radius


def driver():
	[stations, distances, filenames, EQtime, earthquakes_dir, offsets_dir] = configure();
	dataobj_list = inputs(stations, filenames);
	[noeq_objects, east_slope_obj, north_slope_obj] = compute(dataobj_list, distances, earthquakes_dir, offsets_dir, EQtime);
	outputs(noeq_objects, east_slope_obj, north_slope_obj);
	return;


def configure():
	EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
	EQtime  = dt.datetime.strptime("20140310", "%Y%m%d");
	# EQtime  = dt.datetime.strptime("20100110", "%Y%m%d");
	# EQtime  = dt.datetime.strptime("20111017", "%Y%m%d");
	# EQtime  = dt.datetime.strptime("20151017", "%Y%m%d");
	# EQtime  = dt.datetime.strptime("20131017", "%Y%m%d");
	earthquakes_dir = earthquakes_dir="../GPS_POS_DATA/Event_Files/";
	offsets_dir = "../GPS_POS_DATA/Offsets/";
	radius=280;  # km. 
	stations, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius);
	filenames=[];
	for station in stations:
		filenames.append("../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos");
	return [stations, distances, filenames, EQtime, earthquakes_dir, offsets_dir];

def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;


def compute(dataobj_list, distances, earthquakes_dir, offsets_dir, EQtime):

	time_duration = 2; # years

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
		[east_slope_before, north_slope_before, vert_slope_before]=gps_ts_functions.get_slope(newobj,starttime=EQtime-dt.timedelta(days=time_duration*365),endtime=EQtime);
		[east_slope_after, north_slope_after, vert_slope_after]=gps_ts_functions.get_slope(newobj,starttime=EQtime,endtime=EQtime+dt.timedelta(days=time_duration*365));
		east_slope_after=np.round(east_slope_after,decimals=1);
		east_slope_before=np.round(east_slope_before,decimals=1);
		east_slope_obj.append([east_slope_before, east_slope_after]);
		north_slope_after=np.round(north_slope_after,decimals=1);
		north_slope_before=np.round(north_slope_before,decimals=1);
		north_slope_obj.append([north_slope_before, north_slope_after]);

	return [noeq_objects, east_slope_obj, north_slope_obj];


def outputs(noeq_objects, east_slope_obj, north_slope_obj):
	lonW=-125.7
	lonE=-122.0
	latS=39.5
	latN=42.0
	ofile=open('accelerations.txt','w');
	for i in range(len(noeq_objects)):
		ofile.write("%f %f %f %f 0 0 0\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1]-east_slope_obj[i][0], north_slope_obj[i][1]-north_slope_obj[i][0]) );
		# ofile.write("%f %f %f %f 0 0 0 %s\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1]-east_slope_obj[i][0], north_slope_obj[i][1]-north_slope_obj[i][0], noeq_objects[i].name) );
	ofile.close();
	subprocess.call(['./accel_map_gps.gmt',str(lonW),str(lonE),str(latS),str(latN)],shell=False);
	return;

