#!/env/python

"""
Program that takes the difference between two station positions, and makes a time series. 
This shows us how 2 GPS stations are moving relative to one another. 
I have the tenv files downloaded from UNR's website in the data directory. 
"""

import collections
import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt 
from scipy import signal
import haversine as haversine



data_collection = collections.namedtuple("data_collection",['name','coords','yyyymmdd','dtarray','decyeararray','dN', 'dE','dU','Sn','Se','Su']);

def main_program():
	[input_file1, input_file2,coordinates_file] = configure();
	[array1, array2] = inputs(input_file1, input_file2,coordinates_file);
	[matched_dtarray_east, east_diff, matched_dtarray_north, north_diff, matched_dtarray_vert, vert_diff]=compute(array1, array2);
	outputs(array1, array2, matched_dtarray_east, east_diff, matched_dtarray_north, north_diff, matched_dtarray_vert, vert_diff);
	return;


# -------------- CONFIGURE ----------------
def configure():
	data_dir="UNR_Data/"
	station1="CME6"
	station2="P160"
	input_file1=data_dir+station1+".NA12.tenv"
	input_file2=data_dir+station2+".NA12.tenv"
	coordinates_file=data_dir+"coordinates_file.txt"
	return [input_file1, input_file2,coordinates_file];

# -------------- INPUTS ----------------
def inputs(input_file1, input_file2,coordinates_file):
	station_array1=read_UNR_magnet_file(input_file1, coordinates_file);
	station_array2=read_UNR_magnet_file(input_file2, coordinates_file);
	return [station_array1, station_array2];

def read_UNR_magnet_file(filename, coordinates_file):
	[decyeararray,east,north,vert,sig_e,sig_n,sig_v]=np.loadtxt(filename,usecols=(2,8,10,12,14,15,16),skiprows=1,unpack=True);
	yyyymmdd_array=[];
	dtarray=[];
	ifile=open(filename);
	ifile.readline();
	for line in ifile:
		station_name=line.split()[0];
		yyMMMdd=line.split()[1];  # has format 07SEP19
		mydateobject=dt.datetime.strptime(yyMMMdd,"%y%b%d");
		dtarray.append(mydateobject);
		yyyymmdd_array.append(dt.datetime.strftime(mydateobject,"%Y%m%d"));

	coords = get_coordinates_for_station(station_name, coordinates_file);  # format [lon, lat]
	my_data_object=data_collection(name=station_name,coords=coords, yyyymmdd=yyyymmdd_array, dtarray=dtarray, decyeararray=decyeararray, dN=north, dE=east, dU=vert, Sn=sig_n, Se=sig_e, Su=sig_v);
	return my_data_object;

def get_coordinates_for_station(station_name,coordinates_file):
	ifile=open(coordinates_file,'r');
	for line in ifile:
		temp=line.split();
		name=temp[0];
		if name==station_name:
			lon=float(temp[1]);
			lat=float(temp[2]);
			break;
	ifile.close();
	return [lon,lat];



# -------------- COMPUTATION ----------------
def compute(array1, array2):
	
	[matched_dtarray_east, matched_east1, matched_east2]=pair_two_ts(array1.dtarray, array2.dtarray, array1.dE, array2.dE);
	east_diff=np.subtract(matched_east1,matched_east2);
	east_diff=east_diff-np.mean(east_diff);
	
	[matched_dtarray_north, matched_north1, matched_north2]=pair_two_ts(array1.dtarray, array2.dtarray, array1.dN, array2.dN);
	north_diff=np.subtract(matched_north1,matched_north2);
	north_diff=north_diff-np.mean(north_diff);
	
	[matched_dtarray_vert, matched_vert1, matched_vert2]=pair_two_ts(array1.dtarray, array2.dtarray, array1.dU, array2.dU);
	vert_diff=np.subtract(matched_vert1,matched_vert2);	
	vert_diff=vert_diff-np.mean(vert_diff);
	return [matched_dtarray_east, east_diff, matched_dtarray_north, north_diff, matched_dtarray_vert, vert_diff];


def pair_two_ts(dtarray1, dtarray2, data1, data2):  # this function returns the subset of the two time series that are recorded on the same days. 
	matched_dtarray=[];
	matched_data1=[];
	matched_data2=[];
	for i,item in enumerate(dtarray1):
		if item in dtarray2:
			j=dtarray2.index(item);
			matched_dtarray.append(item);
			matched_data1.append(data1[i]);
			matched_data2.append(data2[j]);
	return [matched_dtarray, matched_data1, matched_data2];



# -------------- OUTPUTS ----------------
def outputs(data1, data2, matched_dtarray_east, east_diff, matched_dtarray_north, north_diff, matched_dtarray_vert, vert_diff):

	# Detrending
	east_detrended=signal.detrend(east_diff);
	north_detrended=signal.detrend(north_diff);
	vert_detrended=signal.detrend(vert_diff);

	# The major figure
	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].set_title(data1.name+" / " + data2.name+" relative motion");
	axarr[0].plot_date(matched_dtarray_east, east_diff*1000.0);
	axarr[0].grid('on');
	axarr[0].set_ylabel('east (mm)');
	ax1=axarr[0].twinx();
	ax1.plot_date(matched_dtarray_east, east_detrended*1000,'r');
	ax1.set_ylabel('detrended (mm)')
	

	axarr[1].plot_date(matched_dtarray_north, north_diff*1000.0);
	axarr[1].grid('on');
	axarr[1].set_ylabel('north (mm)');
	ax2=axarr[1].twinx();
	ax2.plot_date(matched_dtarray_north, north_detrended*1000,'r');
	ax2.set_ylabel('detrended (mm)')
	
	axarr[2].plot_date(matched_dtarray_vert, vert_diff*1000.0);	
	axarr[2].grid('on');
	axarr[2].set_ylabel('vertical (mm)')
	ax3=axarr[2].twinx();
	ax3.plot_date(matched_dtarray_vert, vert_detrended*1000,'r');
	ax3.set_ylabel('detrended (mm)')
	plt.savefig(data1.name+"_"+data2.name+"_relative_motion.jpg");

	Distance = haversine.distance([data1.coords[1],data1.coords[0]], [data2.coords[1], data2.coords[0]] );
	print "Distance in km is..."
	print Distance;

	print "Strain Rate from 120mm relative motion in 9 years:"
	print 0.120/22000.0/9.0
	print "strain per year"

	print "Equivalent Gradient in mm/yr/100km:"
	print (120/9.0)*(100/22);


	# Getting the velocities
	[u1,v1]=get_velocities(data1);
	[u2,v2]=get_velocities(data2);

	[ca_lon,ca_lat]=np.loadtxt("california_bdr",unpack=True);
	plt.figure();
	plt.plot(data1.coords[0],data1.coords[1],'.');
	plt.plot(data2.coords[0],data2.coords[1],'.');
	plt.quiver(data1.coords[0], data1.coords[1], u1, v1);
	plt.quiver(data2.coords[0], data2.coords[1], u2, v2);
	plt.plot(ca_lon,ca_lat,'k');
	plt.ylim([39.8,41])
	plt.xlim([-125,-123.5])
	plt.title(data1.name+" "+data2.name+": "+str(Distance)+"km apart")
	plt.savefig(data1.name+"_"+data2.name+"_basic_map.jpg");

	
	# Figure for dataset 1
	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].set_title(data1.name);
	axarr[0].plot_date(data1.dtarray, (data1.dE-np.mean(data1.dE))*1000.0);
	axarr[0].set_xlim([matched_dtarray_east[0],matched_dtarray_east[-1]]);
	axarr[0].grid('on');
	axarr[0].set_ylabel('east (mm)');
	
	axarr[1].plot_date(data1.dtarray, (data1.dN-np.mean(data1.dN))*1000.0);
	axarr[1].set_xlim([matched_dtarray_north[0],matched_dtarray_north[-1]]);
	axarr[1].grid('on');
	axarr[1].set_ylabel('north (mm)');
	
	axarr[2].plot_date(data1.dtarray, (data1.dU-np.mean(data1.dU))*1000.0);
	axarr[2].set_xlim([matched_dtarray_vert[0],matched_dtarray_vert[-1]]);
	axarr[2].grid('on');
	axarr[2].set_ylabel('vertical (mm)')
	plt.savefig(data1.name+"_position.jpg");


	# Figure for dataset 2
	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].set_title(data2.name);
	axarr[0].plot_date(data2.dtarray, (data2.dE-np.mean(data2.dE))*1000.0);
	axarr[0].set_xlim([matched_dtarray_east[0],matched_dtarray_east[-1]]);
	axarr[0].grid('on');
	axarr[0].set_ylabel('east (mm)');
	
	axarr[1].plot_date(data2.dtarray, (data2.dN-np.mean(data2.dN))*1000.0);
	axarr[1].grid('on');
	axarr[1].set_xlim([matched_dtarray_north[0],matched_dtarray_north[-1]]);
	axarr[1].set_ylabel('north (mm)');
	
	axarr[2].plot_date(data2.dtarray, (data2.dU-np.mean(data2.dU))*1000.0);
	axarr[2].grid('on');
	axarr[2].set_xlim([matched_dtarray_vert[0],matched_dtarray_vert[-1]]);
	axarr[2].set_ylabel('vertical (mm)')
	plt.savefig(data2.name+"_position.jpg");

	return;



def get_velocities(data1):
	uparams=np.polyfit(data1.decyeararray, data1.dE,1);
	vparams=np.polyfit(data1.decyeararray, data1.dN,1);
	return [uparams[0], vparams[0]];


