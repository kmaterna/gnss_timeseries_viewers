# Makes a basic python plot of time series position
# Takes in a namedTuple collection with the data. 
# Makes a basic plot. 

import numpy as np
import matplotlib.pyplot as plt 
import collections
import datetime as dt 
import subprocess
from scipy import signal
import gps_io_functions

# For reference of how this gets returned from the read functions.
Timeseries = collections.namedtuple("Timeseries",['name','coords','yyyymmdd','dN', 'dE','dU','Sn','Se','Su']);
Parameters = collections.namedtuple("Parameters",['station','filename','outliers_remove','earthquakes_remove','earthquakes_dir','reference_frame']);

def view_single_station(station_name):
	MyParams=configure(station_name);
	[myData]=gps_io_functions.read_pbo_pos_file(MyParams.filename);
	updatedData=compute(myData,MyParams);
	single_ts_plot(myData);



# -------------- CONFIGURE ------------ # 
def configure(station):
	filename="../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos"
	earthquakes_dir="../GPS_POS_DATA/Event_Files/"
	outliers_remove=0;
	earthquakes_remove=1;
	reference_frame=0;
	MyParams=Parameters(station=station,filename=filename, outliers_remove=outliers_remove, earthquakes_remove=earthquakes_remove, earthquakes_dir=earthquakes_dir, reference_frame=reference_frame);
	return MyParams;



# -------------- COMPUTE ------------ # 
def compute(myData, MyParams):
	if MyParams.earthquakes_remove==1:
		table = subprocess.check_output("grep "+MyParams.station+" "+MyParams.earthquakes_dir+"*kalts.evt",shell=True);
		print(table);
	newData=0;
	return newData;




# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj):
	# Detrending
	east_detrended=signal.detrend(ts_obj.dE);
	north_detrended=signal.detrend(ts_obj.dN);
	vert_detrended=signal.detrend(ts_obj.dU);

	dtarray=[];
	for i in ts_obj.yyyymmdd:
		dtarray.append(dt.datetime.strptime(str(int(i)),"%Y%m%d"));

	# The major figure
	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].set_title(ts_obj.name);
	axarr[0].plot_date(dtarray, ts_obj.dE*1000.0);
	axarr[0].grid('on');
	axarr[0].set_ylabel('east (mm)');
	ax1=axarr[0].twinx();
	ax1.plot_date(dtarray, east_detrended*1000,'r');
	ax1.set_ylabel('detrended (mm)')

	axarr[1].plot_date(dtarray, ts_obj.dN*1000.0);
	axarr[1].grid('on');
	axarr[1].set_ylabel('north (mm)');
	ax2=axarr[1].twinx();
	ax2.plot_date(dtarray, north_detrended*1000,'r');
	ax2.set_ylabel('detrended (mm)')
	
	axarr[2].plot_date(dtarray, ts_obj.dU*1000.0);	
	axarr[2].grid('on');
	axarr[2].set_ylabel('vertical (mm)')
	ax3=axarr[2].twinx();
	ax3.plot_date(dtarray, vert_detrended*1000,'r');
	ax3.set_ylabel('detrended (mm)')
	plt.savefig("single_plots/"+ts_obj.name+"_ts_plot.jpg");

	return;
