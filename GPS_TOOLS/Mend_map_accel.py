# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers and earthquakes. Then identify 
# Step 4: Produce a table and plot of accelerations before/after time ranges. 

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

# Nice Reference: 
# EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
# EQcoords=[-122.834, 37.829]; # San Francisco Bay Area
# EQcoords=[-123.834, 39.029]; # North Bay Area
# EQcoords=[-125.134, 43.829]; # Oregon
# EQcoords=[-125.134, 46.829]; # Washington
# EQcoords=[-125.134, 48.829]; # Canada


import numpy as np 
import datetime as dt 
import subprocess
import gps_io_functions
import gps_ts_functions
import stations_within_radius


def driver(EQcoords, outfile_name, deltat1, deltat2):
	[stations, distances, filenames, earthquakes_dir, offsets_dir, map_coords, dt1_start, dt1_end, dt2_start, dt2_end, outfile_name] = configure(EQcoords, outfile_name, deltat1, deltat2);
	dataobj_list = inputs(stations, filenames);
	[noeq_objects, east_slope_obj, north_slope_obj] = compute(dataobj_list, distances, earthquakes_dir, offsets_dir, dt1_start, dt1_end, dt2_start, dt2_end);
	outputs(noeq_objects, east_slope_obj, north_slope_obj, map_coords,outfile_name);
	return;


def configure(EQcoords, outfile_name, deltat1, deltat2):
	dt1_start  = dt.datetime.strptime(deltat1[0], "%Y%m%d");
	dt1_end  = dt.datetime.strptime(deltat1[1], "%Y%m%d");
	dt2_start  = dt.datetime.strptime(deltat2[0], "%Y%m%d");
	dt2_end  = dt.datetime.strptime(deltat2[1], "%Y%m%d");
	earthquakes_dir = "../../GPS_POS_DATA/PBO_Event_Files/";
	offsets_dir = "../../GPS_POS_DATA/Offsets/";
	radius=450;  # km. 
	map_coords=[EQcoords[0]-0.6, EQcoords[0]+4, EQcoords[1]-2.0, EQcoords[1]+2.0];
	stations, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius, map_coords);
	filenames=[];
	for station in stations:
		filenames.append("../../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos");
	return [stations, distances, filenames, earthquakes_dir, offsets_dir, map_coords, dt1_start, dt1_end, dt2_start, dt2_end, outfile_name];


def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;


def compute(dataobj_list, distances, earthquakes_dir, offsets_dir, dt1_start, dt1_end, dt2_start, dt2_end):
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
		newobj=gps_ts_functions.make_detrended_option(newobj, 1, 'fit');
		noeq_objects.append(newobj);

		# Get the pre-event and post-event velocities (earthquakes removed)
		[east_slope_before, north_slope_before, vert_slope_before]=gps_ts_functions.get_slope(newobj,starttime=dt1_start,endtime=dt1_end);
		[east_slope_after, north_slope_after, vert_slope_after]=gps_ts_functions.get_slope(newobj,starttime=dt2_start,endtime=dt2_end);
		if east_slope_before+east_slope_after+north_slope_before+north_slope_after+vert_slope_before+vert_slope_after==np.nan:
			east_slobe_obj.append(0);
			north_slope_obj.append(0);
		else:
			east_slope_after=np.round(east_slope_after,decimals=1);
			east_slope_before=np.round(east_slope_before,decimals=1);
			east_slope_obj.append([east_slope_before, east_slope_after]);
			north_slope_after=np.round(north_slope_after,decimals=1);
			north_slope_before=np.round(north_slope_before,decimals=1);
			north_slope_obj.append([north_slope_before, north_slope_after]);
	return [noeq_objects, east_slope_obj, north_slope_obj];


def outputs(noeq_objects, east_slope_obj, north_slope_obj, map_coords,outfile_name):
	# basename=outfile_name.split(".")[0]
	ofile=open('accelerations.txt','w');
	for i in range(len(noeq_objects)):
		ofile.write("%f %f %f %f 0 0 0\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1]-east_slope_obj[i][0], north_slope_obj[i][1]-north_slope_obj[i][0]) );
		# ofile.write("%f %f %f %f 0 0 0 %s\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1]-east_slope_obj[i][0], north_slope_obj[i][1]-north_slope_obj[i][0], noeq_objects[i].name) );
	ofile.close();
	subprocess.call(['./accel_map_gps.gmt',str(map_coords[0]),str(map_coords[1]),str(map_coords[2]),str(map_coords[3]),outfile_name],shell=False);
	print('./accel_map_gps.gmt '+str(map_coords[0])+' '+str(map_coords[1])+' '+str(map_coords[2])+' '+str(map_coords[3])+' '+outfile_name);
	return;

