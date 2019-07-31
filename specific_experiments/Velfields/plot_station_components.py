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
	network = 'pbo';  # choices: unr, pbo, cwu, nmt, nldas, gldas, lsdm
	refname = 'NA'; # choices: NA, ITRF
	seasonal_types = ['none','lssq','notch','nldas','lsdm','gldas','grace']; # OPTIONS: none, lssq, notch, stl, grace, gldas, nldas, lsdm, shasta, oroville. 
	station = "P170"
	return [station, network, refname, seasonal_types];

def inputs(station_name, network, refframe):
	[dataobj, offsetobj, eqobj] = gps_input_pipeline.get_station_data(station_name, network, refframe);
	return [dataobj, offsetobj, eqobj];

def compute(dataobj, offsetobj, eqobj, seasonal_types):
	trended_objects=[];

	# Remove the earthquakes
	newobj=offsets.remove_offsets(dataobj,offsetobj);
	newobj=offsets.remove_offsets(newobj, eqobj);

	for i in seasonal_types:

		if i == 'none':
			trended=newobj;
		else:
			if newobj.name=='P349':
				newobj_detrended, newobj=gps_seasonal_removals.make_detrended_ts(newobj,1,'shasta', remove_trend=0);  # special considerations for stations near lakes
			if newobj.name=='ORVB':
				newobj_detrended, newobj=gps_seasonal_removals.make_detrended_ts(newobj,1,'oroville', remove_trend=0);
			# For each object, remove seasonals in the specified way. 
			trended=gps_seasonal_removals.make_detrended_ts(newobj, 1, i, remove_trend=0);  # remove seasonals but do not remove the trend. 

		trended_objects.append(trended);

	return [trended_objects];


def write_outputs(trended_objects, seasonal_types):

	colors=['black','blue','red','green','purple','orange','cyan'];
	offset=20; 

	plt.figure(figsize=(16,3));
	for i in range(len(seasonal_types)):
		plot_function = [offset*i + trended_objects[i].dE[j] for j in range(len(trended_objects[i].dE))];
		plt.plot(trended_objects[i].dtarray, plot_function,'.',color=colors[i],label=seasonal_types[i]);
	plt.legend(loc=3);
	plt.savefig(trended_objects[0].name+"_east_no_seasonals");

	plt.figure(figsize=(15,15));
	for i in range(len(seasonal_types)):
		plot_function = [offset*i + trended_objects[i].dU[j] for j in range(len(trended_objects[i].dU))];
		plt.plot(trended_objects[i].dtarray, plot_function,'.',color=colors[i],label=seasonal_types[i]);
	plt.ylim([-20, 180])
	plt.legend(loc=3,fontsize=16);
	plt.savefig(trended_objects[0].name+"_vert_no_seasonals");

	return;


if __name__=="__main__":
	[station_name, network, refframe, seasonal_types] = configure();
	[dataobj, offsetobj, eqobj] = inputs(station_name, network, refframe);
	[trended_objects] = compute(dataobj, offsetobj, eqobj, seasonal_types);
	write_outputs(trended_objects, seasonal_types);

