# Makes a basic python plot of time series position
# Takes in a namedTuple collection with the data. 
# Makes a basic plot. 

import numpy as np
import matplotlib.pyplot as plt 
import collections
import datetime as dt 
import subprocess, sys
from scipy import signal
import gps_io_functions

# For reference of how this gets returned from the read functions.
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Parameters = collections.namedtuple("Parameters",['station','filename','outliers_remove','earthquakes_remove','earthquakes_dir','reference_frame']);

def view_single_station(station_name):
	MyParams=configure(station_name);
	[myData]=gps_io_functions.read_pbo_pos_file(MyParams.filename);
	updatedData=compute(myData,MyParams);
	single_ts_plot(updatedData,MyParams);



# -------------- CONFIGURE ------------ # 
def configure(station):
	filename="../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos"
	earthquakes_dir="../GPS_POS_DATA/Event_Files/"
	earthquakes_remove = 1;
	outliers_remove    = 0;	
	outliers_def       = 15.0;  # mm away from average. 
	reference_frame    = 0;
	MyParams=Parameters(station=station,filename=filename, outliers_remove=outliers_remove, earthquakes_remove=earthquakes_remove, earthquakes_dir=earthquakes_dir, reference_frame=reference_frame);
	return MyParams;



# -------------- COMPUTE ------------ # 
def compute(myData, MyParams):
	# First step: remove earthquakes
	newData=myData; 

	if MyParams.earthquakes_remove==1:
		newData=remove_earthquakes(myData, MyParams);

	# Second step: remove outliers
	if MyParams.outliers_remove==1:
		newData=remove_outliers(newData, MyParams);

	return newData;




def remove_earthquakes(Data0, MyParams):
	# Building earthquake table
	table = subprocess.check_output("grep "+MyParams.station+" "+MyParams.earthquakes_dir+"*kalts.evt",shell=True);
	print("building earthquake table...")
	print(table);
	e_offset=[]; n_offset=[]; u_offset=[]; evdt=[];
	tablesplit=table.split('\n');
	for item in tablesplit:  # for each earthquake
		if len(item)==0:
			continue;  # if we're at the end, move on. 
		words = item.split();
		filename=words[0];
		e_offset.append(float(words[3]));  # in m
		n_offset.append(float(words[4]));
		u_offset.append(float(words[8]));
		evdate = filename.split('/')[-1];
		evdate = evdate[4:10];
		year=evdate[0:2];
		month=evdate[2:4];
		day=evdate[4:6];
		year="20"+year;
		evdt.append(dt.datetime.strptime(year+month+day,"%Y%m%d"));

	newdtarray=[]; newdN=[]; newdE=[]; newdU=[];
	# Removing offsets
	for i in range(len(Data0.dtarray)):
		# For each day...
		tempE=Data0.dE[i];
		tempN=Data0.dN[i];
		tempU=Data0.dU[i];
		for j in range(len(evdt)):
			# print("removing %f mm from east at %s" % (e_offset[j], evdt[j]));
			if Data0.dtarray[i]>=evdt[j]:
				tempE=tempE-e_offset[j];
				tempN=tempN-n_offset[j];
				tempU=tempU-u_offset[j];
		newdtarray.append(Data0.dtarray[i]);
		newdE.append(tempE);
		newdN.append(tempN);
		newdU.append(tempU);
	
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=evdt);
	return newData;


def remove_outliers(myData, MyParams):
	newData=[];
	return newData;




# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, MyParams):
	# Detrending
	east_detrended=signal.detrend(ts_obj.dE);
	north_detrended=signal.detrend(ts_obj.dN);
	vert_detrended=signal.detrend(ts_obj.dU);

	# The major figure
	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].set_title(ts_obj.name);
	axarr[0].plot_date(ts_obj.dtarray, ts_obj.dE);
	axarr[0].grid('on');
	axarr[0].set_ylabel('east (mm)');
	bottom,top=axarr[0].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[0].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k');
	ax1=axarr[0].twinx();
	ax1.plot_date(ts_obj.dtarray, east_detrended,'r');
	ax1.set_ylabel('detrended (mm)')

	axarr[1].plot_date(ts_obj.dtarray, ts_obj.dN);
	axarr[1].grid('on');
	axarr[1].set_ylabel('north (mm)');
	bottom,top=axarr[1].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[1].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k');	
	ax2=axarr[1].twinx();
	ax2.plot_date(ts_obj.dtarray, north_detrended,'r');
	ax2.set_ylabel('detrended (mm)')
	
	axarr[2].plot_date(ts_obj.dtarray, ts_obj.dU);	
	axarr[2].grid('on');
	axarr[2].set_ylabel('vertical (mm)')
	bottom,top=axarr[2].get_ylim();
	for i in range(len(ts_obj.EQtimes)):
		axarr[2].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k');	
	ax3=axarr[2].twinx();
	ax3.plot_date(ts_obj.dtarray, vert_detrended,'r');
	ax3.set_ylabel('detrended (mm)')
	

	if MyParams.earthquakes_remove==1:
		plt.savefig("single_plots/"+ts_obj.name+"_ts_noeq.jpg");
	else:
		plt.savefig("single_plots/"+ts_obj.name+"_ts.jpg");

	return;
