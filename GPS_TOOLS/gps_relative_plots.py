# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and steps. Make the stations relative to a base station.
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
import outputs_gps_stacks
import remove_ets_events


Parameters=collections.namedtuple("Parameters",['expname','proc_center','refframe','center','radius','stations','distances','blacklist','ref_station','outdir', 'outname']);


def driver():
	myparams = configure();
	[dataobj_list, offsetobj_list, eqobj_list, paired_distances] = gps_input_pipeline.multi_station_inputs(myparams.stations, myparams.blacklist, myparams.proc_center, myparams.refframe, myparams.distances);

	[relative_objects, sorted_distances] = compute(dataobj_list, offsetobj_list, eqobj_list, paired_distances, myparams.ref_station);

	outputs_gps_stacks.horizontal_full_ts(relative_objects, sorted_distances, myparams, "noeq_noseasons");
	outputs_gps_stacks.vertical_full_ts(relative_objects, sorted_distances, myparams);
	outputs_gps_stacks.horizontal_filtered_plots(relative_objects, sorted_distances, myparams);
	outputs_gps_stacks.vertical_filtered_plots(relative_objects, sorted_distances, myparams);
	return;


def configure():
	# center=[-119.0, 37.7];     expname='LVC';  radius = 30; # km
	# center=[-115.5, 32.85]; expname='SSGF'; radius = 20; 
	# center=[-115.5, 33]; expname='SSGF'; radius = 15; 
	center=[-115.5, 33]; expname='SSGF'; radius =40; 

	proc_center='nmt';   # WHICH DATASTREAM DO YOU WANT?
	refframe = 'NA';     # WHICH REFERENCE FRAME? 

	stations, distances = stations_within_radius.get_stations_within_radius(center, radius, network=proc_center);
	blacklist=["P316","P170","P158","TRND","P203","BBDM","KBRC","RYAN","BEAT","CAEC","MEXI","BOMG","FSHB"];  # This is global, just keeps growing
	reference_station = "P502"
	outdir="Rel_"+reference_station+"_"+expname+"_"+proc_center
	subprocess.call(["mkdir","-p",outdir],shell=False);
	outname=expname+"_"+str(center[0])+"_"+str(center[1])+"_"+str(radius)
	myparams=Parameters(expname=expname, proc_center=proc_center, refframe=refframe, center=center, radius=radius, stations=stations, distances=distances, blacklist=blacklist, ref_station=reference_station, outdir=outdir, outname=outname);
	return myparams;



def compute(dataobj_list, offsetobj_list, eqobj_list, distances, ref_station):

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_offsets = [x for _,x in sorted(zip(latitudes_list, offsetobj_list))];  # the raw, sorted data. 
	sorted_eqs = [x for _,x in sorted(zip(latitudes_list, eqobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.


	no_offset_objects = [];
	relative_objects = [];

	# Detrended objects (or objects with trends and no offsets; depends on what you want.)
	for i in range(len(sorted_objects)):
		newobj=offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
		newobj=offsets.remove_offsets(newobj,sorted_eqs[i]);
		no_offset_objects.append(newobj);  # still has trends, doesn't have offsets
		if newobj.name==ref_station:
			ref_index=i;

	for i in range(len(sorted_objects)):
		newobj=gps_ts_functions.get_referenced_data(no_offset_objects[i], no_offset_objects[ref_index]);
		relative_objects.append(newobj)


	return [relative_objects, sorted_distances];


