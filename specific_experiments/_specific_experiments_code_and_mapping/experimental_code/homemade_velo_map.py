# Python viewing of velocities

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and offsets
# Step 4: Plot maps of vertical velocity

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
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
import pygmt


Parameters=collections.namedtuple("Parameters",['expname','proc_center','refframe','center','radius','stations','distances','blacklist','outdir', 'outname']);


def driver():
	myparams = configure();
	[dataobj_list, offsetobj_list, eqobj_list, temp] = gps_input_pipeline.multi_station_inputs(myparams.stations, myparams.blacklist, myparams.proc_center, myparams.refframe, myparams.distances);
	[east_slope_list, north_slope_list, vert_slope_list] = compute(dataobj_list, offsetobj_list, eqobj_list);
	pygmt_vertical_map(myparams, dataobj_list, vert_slope_list);
	return;

def configure():
	# center=[-124, 40.5]; expname='Humboldt'; radius = 100; # 
	# center=[-122.5, 40.5]; expname='Chico'; radius = 75; # 
	# center=[-119.0, 37.7];     expname='LVC';  radius = 30; # km
	center=[-115.5, 33]; expname='SSGF'; radius = 80; 

	proc_center='unr';   # WHICH DATASTREAM DO YOU WANT?
	refframe = 'NA';     # WHICH REFERENCE FRAME? 

	stations, lons, lats, distances = stations_within_radius.get_stations_within_radius(center, radius, network=proc_center);
	blacklist=["P340","P316","P170","P158","TRND","P203","BBDM","KBRC","RYAN","BEAT","CAEC","MEXI"];  # This is global, just keeps growing
	outdir=expname+"_"+proc_center
	subprocess.call(["mkdir","-p",outdir],shell=False);
	outname=expname+"_"+str(center[0])+"_"+str(center[1])+"_"+str(radius)
	myparams=Parameters(expname=expname, proc_center=proc_center, refframe=refframe, center=center, radius=radius, stations=stations, distances=distances, blacklist=blacklist, outdir=outdir, outname=outname);
	return myparams;


def compute(dataobj_list, offsetobj_list, eqobj_list):
	east_slope_list=[]; north_slope_list=[]; vert_slope_list=[];

	# Objects with no earthquakes or seasonals
	for i in range(len(dataobj_list)):

		# Remove the steps earthquakes
		newobj=offsets.remove_offsets(dataobj_list[i], offsetobj_list[i]);
		newobj=offsets.remove_offsets(newobj,eqobj_list[i]);
		[east_slope, north_slope, vert_slope, east_std, north_std, vert_std] = gps_ts_functions.get_slope(newobj);
		east_slope_list.append(east_slope)
		north_slope_list.append(north_slope)
		vert_slope_list.append(vert_slope)

	return [east_slope_list, north_slope_list, vert_slope_list];


def pygmt_vertical_map(myparams, ts_objects, vert_slopes):

	offset=0.2;
	geothermals_x=[-115.528300, -115.250000, -115.515300, -115.600000];
	geothermals_y=[32.716700, 32.783300, 33.015300, 33.200000];

	lons=[]; lats=[]; names=[];
	for i in range(len(ts_objects)):
		lons.append(ts_objects[i].coords[0]);
		lats.append(ts_objects[i].coords[1]);
		names.append(ts_objects[i].name);
	region=[min(lons)-offset,max(lons)+offset,min(lats)-offset,max(lats)+offset];		

	# Make a new color bar
	min_vert=np.min(vert_slopes);
	if min_vert < -8:
		min_vert=-8;
	max_vert=np.max(vert_slopes);
	if max_vert-min_vert < 5.0:
		label_interval=0.5;
	else:
		label_interval=2.0;
	pygmt.makecpt(C="jet",T=str(min_vert-0.1)+"/"+str(max_vert+0.1)+"/0.1",H="mycpt.cpt");

	fig = pygmt.Figure()
	fig.basemap(region=region,projection="M8i",B="0.25");
	fig.coast(shorelines="0.5p,black",G='peachpuff2',S='skyblue',D="h");
	fig.coast(N='1',W='1.0p,black');
	fig.coast(N='2',W='0.5p,black');
	fig.text(x=[i+0.035 for i in lons],y=lats,text=names,font='15p,Helvetica-Bold,black');
	fig.plot(x=geothermals_x,y=geothermals_y,S='i0.2i',G="purple",W='0.5p,black');
	fig.plot(x=lons,y=lats,S='c0.2i',C="mycpt.cpt",G=vert_slopes,W='0.5p,black');
	fig.plot(x=myparams.center[0],y=myparams.center[1],S='a0.1i',G='red',W='0.5p,red')
	fig.colorbar(D="JBC+w4.0i+h",C="mycpt.cpt",G=str(min_vert)+"/"+str(max_vert),B=["x"+str(label_interval),"y+LVelocity(mm/yr)"])
	fig.savefig(myparams.outdir+"/"+myparams.outname+'_map.png');
	return;


