# Velocity field experiment
# July 2019
# Make velocity fields from various GPS solutions. 
# Removing earthquakes, seasonals, etc. 
# How different are they really? 
# 2D vs 3D? 



import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import glob
import subprocess, sys
import gps_io_functions
import gps_ts_functions
import gps_seasonal_removals
import gps_input_pipeline
import offsets
import stations_within_radius
import haversine
import grace_ts_functions
import remove_ets_events



def configure():
	overall_size = 'medium';
	network = 'pbo';  # choices: unr, pbo, cwu, nmt, nldas, gldas, lsdm
	refname = 'NA'; # choices: NA, ITRF
	seasonal_type = 'notch'; # OPTIONS: none, lssq, notch, stl, grace, gldas, nldas, lsdm, shasta, oroville. 
	outfile = "Fields/"+network+"_"+refname+"_"+seasonal_type+"_"+"small"+"_velocities.txt";

	if overall_size=='small':
		map_coords=[-125.6, -123.0, 39.5, 41.5];
	elif overall_size=='medium':
		map_coords=[-125.6, -120.0, 38.0, 43.0];
	elif overall_size=='huge':
		map_coords=[-125.6, -110.0, 32.5, 48.5];
	else:
		map_coords=[-125.6, -121.0, 39.0, 41.5];
	stations = stations_within_radius.get_stations_within_box(map_coords, network);
	print("Box is");
	print(map_coords);
	print("Returning %d stations and putting their velocities in %s " % (len(stations), outfile));
	return [stations, network, refname, seasonal_type, outfile];

def inputs(station_names, network, refframe):
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[];
	for station_name in station_names:
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, network, refframe);
		if myData==[]:
			continue;
		else:
			dataobj_list.append(myData);
			offsetobj_list.append(offset_obj);
			eqobj_list.append(eq_obj);
	return [dataobj_list, offsetobj_list, eqobj_list];

def compute(dataobj_list, offsetobj_list, eqobj_list, seasonal_type):

	noeq_objects = [];	# No earthquakes objects
	east_slope_obj=[]; north_slope_obj=[]; vert_slope_obj=[]; 

	# The main processing loop for slopes. 
	for i in range(len(dataobj_list)):
		# Remove the earthquakes
		print(dataobj_list[i].name);
		newobj=offsets.remove_offsets(dataobj_list[i],offsetobj_list[i]);
		newobj=offsets.remove_offsets(newobj, eqobj_list[i]);
		if seasonal_type == 'none':
			newobj=newobj;
		else:
			if newobj.name=='P349':
				newobj=gps_seasonal_removals.make_detrended_ts(newobj,1,'shasta', remove_trend=0);  # special considerations for stations near lakes
			if newobj.name=='ORVB':
				newobj=gps_seasonal_removals.make_detrended_ts(newobj,1,'oroville', remove_trend=0);
			# For each object, remove seasonals in the specified way. 
			newobj=gps_seasonal_removals.make_detrended_ts(newobj, 1, seasonal_type, remove_trend=0);  # remove seasonals but do not remove the trend. 

		noeq_objects.append(newobj);
		
		# Get slopes
		[east_slope, north_slope, vert_slope, esig0, nsig0, usig0]=gps_ts_functions.get_slope(newobj);
		print(east_slope);
		east_slope_obj.append([east_slope, esig0]);
		north_slope_obj.append([north_slope, nsig0]);
		vert_slope_obj.append([vert_slope, usig0]);

	return [noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj]; 

def write_outputs(noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj, outfile):
	ofile=open(outfile,'w');
	for i in range(len(noeq_objects)):
		ofile.write("%f %f %f %f %f %f %f %f 0 %s\n" %(noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][0], north_slope_obj[i][0], vert_slope_obj[i][0], east_slope_obj[i][1], north_slope_obj[i][1], vert_slope_obj[i][1], noeq_objects[i].name) );
	ofile.close();
	return;


if __name__=="__main__":
	[station_names, network, refframe, seasonal_type, outfile] = configure();
	[dataobj_list, offsetobj_list, eqobj_list] = inputs(station_names, network, refframe);
	[noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj] = compute(dataobj_list, offsetobj_list, eqobj_list, seasonal_type); 
	write_outputs(noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj, outfile);

