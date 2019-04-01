#!/usr/bin/env python

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import collections
import gps_input_pipeline
import offsets
import gps_ts_functions
import gps_seasonal_removals

Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


def configure():
	stations=["P349","P060","P348","P338","P341","WDCB","ORVB"];
	lakes = ["shasta","shasta","shasta","shasta","shasta","shasta","oroville"];
	time0=dt.datetime.strptime("20100610","%Y%m%d");
	time1=dt.datetime.strptime("20140317","%Y%m%d");
	time2=dt.datetime.strptime("20161208","%Y%m%d");
	time3=dt.datetime.strptime("20180915","%Y%m%d");
	outfile='OUTPUTS/model_vs_data.txt'
	return [stations, lakes, time0, time1, time2, time3, outfile];

def inputs(stations, lakes):
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[];
	for station_name in stations:
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, 'pbo');
		dataobj_list.append(myData);
		offsetobj_list.append(offset_obj);
		eqobj_list.append(eq_obj);
	[loadingobj_list] = input_loading_manystations(stations, lakes);
	return [dataobj_list, offsetobj_list, eqobj_list, loadingobj_list];

def input_loading_manystations(stations, lakes):
	loadingobj_list=[];
	for i in range(len(stations)):
		loading_obj=input_loading_onestation(stations[i], lakes[i]);
		loadingobj_list.append(loading_obj);
	return [loadingobj_list];

def input_loading_onestation(station, lake):
	ifile=open("../../GPS_POS_DATA/Lake_loading/"+station+"_"+lake+"_defo.txt",'r');
	dtarray = []; u = []; v = []; w = []; 
	for line in ifile:
		temp=line.split();
		dtarray.append(dt.datetime.strptime(temp[0],"%Y-%m-%d"));
		u.append(float(temp[4]));
		v.append(float(temp[5]));
		w.append(float(temp[6]));
	# HERE WE WILL MAKE A NEW DATA OBJECT
	ifile.close();
	loading_defo = Timeseries(name=station, coords=[], dtarray=dtarray, dE=u, dN=v, dU=w, Sn=[], Se=[], Su=[], EQtimes=[]);
	return loading_defo;

def compute(dataobj_list, offsetobj_list, eqobj_list, loadingobj_list, time0, time1, time2, time3):
	# Initialize output objects
	noeq_objects = []; 
	for i in range(len(dataobj_list)):
		# Remove the earthquakes and offsets
		newobj=offsets.remove_offsets(dataobj_list[i], offsetobj_list[i]);
		newobj=offsets.remove_offsets(newobj,eqobj_list[i]);
		newobj=gps_ts_functions.remove_outliers(newobj,15);  # 15mm horizontal outliers
		newobj=gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq');  # can remove seasonals a few ways
		noeq_objects.append(newobj);

	T3_T2_gps_slopes=[]; T4_T3_gps_slopes=[];
	for i in range(len(noeq_objects)):
		# T2, T3, T4
		[east_T2, north_T2, _, east_std_T2, north_std_T2, _] = gps_ts_functions.get_slope(noeq_objects[i], starttime=time0, endtime=time1);
		[east_T3, north_T3, _, east_std_T3, north_std_T3, _] = gps_ts_functions.get_slope(noeq_objects[i], starttime=time1, endtime=time2);
		[east_T4, north_T4, _, east_std_T4, north_std_T4, _] = gps_ts_functions.get_slope(noeq_objects[i], starttime=time2, endtime=time3);
		T3_T2_gps_slopes.append([east_T3-east_T2, north_T3-north_T2]);
		T4_T3_gps_slopes.append([east_T4-east_T3, north_T4-north_T3]);


	T3_T2_model_slopes=[]; T4_T3_model_slopes=[];
	for i in range(len(loadingobj_list)):
		[east_T2, north_T2, _, east_std_T2, north_std_T2, _] = gps_ts_functions.get_slope(loadingobj_list[i], starttime=time0, endtime=time1);
		[east_T3, north_T3, _, east_std_T3, north_std_T3, _] = gps_ts_functions.get_slope(loadingobj_list[i], starttime=time1, endtime=time2);
		[east_T4, north_T4, _, east_std_T4, north_std_T4, _] = gps_ts_functions.get_slope(loadingobj_list[i], starttime=time2, endtime=time3);
		T3_T2_model_slopes.append([east_T3-east_T2, north_T3-north_T2]);
		T4_T3_model_slopes.append([east_T4-east_T3, north_T4-north_T3]);		

	return [T3_T2_gps_slopes, T4_T3_gps_slopes, T3_T2_model_slopes, T4_T3_model_slopes];

def outputs(T3_T2_gps_slopes, T4_T3_gps_slopes, T3_T2_model_slopes, T4_T3_model_slopes, gps_data, outfile):
	ofile=open(outfile,'w');
	ofile.write("# lon, lat, T3-T2_gps(E,N) T3-T2_model(E,N) T4-T3_gps(E,N) T4-T3_model(E,N)\n");
	for i in range(len(gps_data)):
		ofile.write( "%f %f %f %f %f %f %f %f %f %f\n" % (gps_data[i].coords[0], gps_data[i].coords[1], T3_T2_gps_slopes[i][0], T3_T2_gps_slopes[i][1], T3_T2_model_slopes[i][0], T3_T2_model_slopes[i][1], T4_T3_gps_slopes[i][0], T4_T3_gps_slopes[i][1], T4_T3_model_slopes[i][0], T4_T3_model_slopes[i][1] ) );
	ofile.close();
	return;


if __name__=="__main__":
	[stations, lakes, time0, time1, time2, time3, outfile] = configure();
	[myData, offset_obj, eq_obj, loadingobj_list] = inputs(stations, lakes);
	[T3_T2_gps_slopes, T4_T3_gps_slopes, T3_T2_model_slopes, T4_T3_model_slopes] = compute(myData, offset_obj, eq_obj, loadingobj_list, time0, time1, time2, time3);
	outputs(T3_T2_gps_slopes, T4_T3_gps_slopes, T3_T2_model_slopes, T4_T3_model_slopes, myData, outfile);
