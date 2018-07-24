# Python script to find inflection points in noisy time series data. 
# July 24, 2018
# Kathryn Materna

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt 
from scipy.signal import butter, filtfilt
import subprocess, sys
import gps_io_functions
import gps_ts_functions
import stations_within_radius


# Referecnce : 
# Timeseries = collections.namedtuple("Timeseries",
#      ['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);


def driver(eqtime):
	[stations, filenames, earthquakes_dir, offsets_dir, start_time, end_time, map_coords, outfile_dir] = configure(eqtime);
	dataobj_list = inputs(stations, filenames);
	[noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time]=compute(dataobj_list, earthquakes_dir, offsets_dir, start_time, end_time);
	outputs(noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, outfile_dir);
	return;


# ------------ CONFIGURE ------------- # 
def configure(eqtime):
	EQtime  = dt.datetime.strptime(eqtime, "%Y%m%d");
	pbo_velfile="../../GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt";
	earthquakes_dir = "../../GPS_POS_DATA/PBO_Event_Files/";
	offsets_dir = "../../GPS_POS_DATA/Offsets/";
	pre_event_duration = 3; # years
	post_event_duration = 3; # years
	start_time=EQtime-dt.timedelta(days=pre_event_duration*365);
	end_time=EQtime+dt.timedelta(days=post_event_duration*365);
	map_coords=[-125, -122, 39, 41.5];
	stations = stations_within_radius.get_stations_within_box(map_coords,pbo_velfile);
	filenames=[];
	for station in stations:
		filenames.append("../../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos");
	outfile_dir='Outputs/'
	return [stations, filenames, earthquakes_dir, offsets_dir, start_time, end_time, map_coords, outfile_dir];


# ------------ INPUTS  ------------- # 
def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
		break;  # just for debugging
	return dataobj_list;


# ------------ COMPUTE ------------- # 
def compute(dataobj_list, earthquakes_dir, offsets_dir, start_time, end_time):
	
	# Initialize output objects
	noeq_objects = []; short_objects = [];
	east_filt=[]; north_filt=[]; vert_filt=[]; east_inf_time=[]; north_inf_time=[]; vert_inf_time=[];

	for i in range(len(dataobj_list)):
		# Remove the earthquakes
		newobj=gps_ts_functions.remove_offsets(dataobj_list[i],offsets_dir);
		newobj=gps_ts_functions.remove_earthquakes(newobj,earthquakes_dir);
		newobj=gps_ts_functions.remove_outliers(newobj,15);  # 15mm horizontal outliers
		newobj=gps_ts_functions.detrend_data(newobj);
		newobj=gps_ts_functions.remove_annual_semiannual(newobj);
		short_obj = gps_ts_functions.impose_time_limits(newobj,start_time,end_time);
		# short_obj = newobj;
		noeq_objects.append(newobj);
		short_objects.append(short_obj);

		# Get the inflection points in the timeseries
		float_dtarray = gps_ts_functions.get_float_times(short_obj.dtarray);
		[east_filtered, e_inflection_time]=inflection_with_butterworth(float_dtarray, short_obj.dE);
		[north_filtered, n_inflection_time]=inflection_with_butterworth(float_dtarray, short_obj.dN);
		[vert_filtered, v_inflection_time]=inflection_with_butterworth(float_dtarray, short_obj.dU);

		east_filt.append(east_filtered);
		north_filt.append(north_filtered);
		vert_filt.append(vert_filtered);
		east_inf_time.append(e_inflection_time);
		north_inf_time.append(n_inflection_time);
		vert_inf_time.append(v_inflection_time);

	return [noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time];


# The butterworth filter. 
def inflection_with_butterworth(x, y):

	# Butterworth parameters
	N=3;
	Wn=1/365.0;  # 1/period of cutoff frequency. 

	# East Component
	[b,a]=butter(N, Wn, btype='low',output='ba');
	y_filtered = filtfilt(b, a, y);

	return [y_filtered, 0];



# ------------ OUTPUTS ------------- # 
def outputs(noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, outfile_dir):
	ofile=open(outfile_dir+'inflections.txt','w');
	for i in range(len(noeq_objects)):
		ofile.write("%f %f %f %f %f\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_inf_time[i], north_inf_time[i], vert_inf_time[i]) );
	ofile.close();
	for i in range(len(noeq_objects)):
		output_plots(noeq_objects[i], short_objects[i], east_filt[i], north_filt[i], vert_filt[i], east_inf_time[i], north_inf_time[i], vert_inf_time[i], outfile_dir);
	return;

def output_plots(noeq_obj, short_obj, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, outfile_dir):

	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].plot_date(noeq_obj.dtarray,noeq_obj.dE);
	axarr[0].plot_date(short_obj.dtarray,short_obj.dE,'r.');
	axarr[0].plot_date(short_obj.dtarray,east_filt,'g',linewidth=3);
	axarr[0].set_title(noeq_obj.name+' East');
	
	axarr[1].plot_date(noeq_obj.dtarray,noeq_obj.dN);
	axarr[1].plot_date(short_obj.dtarray,short_obj.dN,'r.');
	axarr[1].plot_date(short_obj.dtarray,north_filt,'g',linewidth=3);
	axarr[1].set_title(noeq_obj.name+' North');
	
	axarr[2].plot_date(noeq_obj.dtarray,noeq_obj.dU);
	axarr[2].plot_date(short_obj.dtarray,short_obj.dU,'r.');
	axarr[2].plot_date(short_obj.dtarray,vert_filt,'g',linewidth=3);
	axarr[2].set_title(noeq_obj.name+' Vertical');

	plt.savefig('test_plot.png');

	return;



# --------- DRIVER ---------- # 
if __name__=="__main__":
	
	eqtime="20140310"; # 2014 M6.8 Earthquake
	# eqtime="20161208"; # # 2016 M6.6 Earthquake
	# eqtime="20100110"; # # 2010 M6.5 Earthquake
	driver(eqtime);

