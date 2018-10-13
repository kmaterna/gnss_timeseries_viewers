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
import gps_seasonal_removals
import gps_input_pipeline
import offsets
import stations_within_radius


# Reference : 
# Timeseries = collections.namedtuple("Timeseries",
#      ['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);


def driver(eqtime, starttime, endtime):
	[stations, fit_table, grace_dir,start_time_infl, end_time_infl,start_time_velo, end_time_velo, N, Wn, seasonal_type, map_coords, outfile_dir] = configure(eqtime, starttime, endtime);
	[dataobj_list, offsetobj_list, eqobj_list] = inputs(stations);
	[noeq_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, east_change, north_change, vert_change]=compute(dataobj_list, offsetobj_list, eqobj_list, fit_table, grace_dir,start_time_infl, end_time_infl,start_time_velo, end_time_velo, seasonal_type, N, Wn);
	outputs(noeq_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, east_change, north_change, vert_change, start_time_infl, end_time_infl, outfile_dir);
	return;


# ------------ CONFIGURE ------------- # 
def configure(eqtime,starttime,endtime):
	EQtime  = dt.datetime.strptime(eqtime, "%Y%m%d");
	fit_table="../../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt"
	grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"
	seasonal_type="lssq";

	if eqtime[0:4]=='2016' or eqtime[0:4]=='2017':
		pre_event_duration = 0.5; # years
		post_event_duration = 0.5; # years
		print('searching %f years, centered on the event time' % ( pre_event_duration+post_event_duration ) );
		# Presently, I'm searching 1 year surrounding the event time for an inflection point. 	
	else:
		pre_event_duration = 0.5; # years
		post_event_duration = 0.5; # years
		print('searching %f years, centered on the event time' % ( pre_event_duration+post_event_duration) );

	start_time_infl=EQtime-dt.timedelta(days=pre_event_duration*365);
	end_time_infl=EQtime+dt.timedelta(days=post_event_duration*365);
	start_time_velo=dt.datetime.strptime(starttime, "%Y%m%d");
	end_time_velo=dt.datetime.strptime(endtime, "%Y%m%d");

	# Butterworth parameters
	N=3;  # Order of butterworth filter
	Wn=1/365.0;  # 1/period (days) of cutoff frequency. 	

	# map_coords=[-125, -118, 36.5, 42.0];  # northern CA
	map_coords=[-125, -110, 32.5, 48.5]; # western US
	# map_coords=[-122, -119, 37.0, 38.0]; # MTJ
	stations = stations_within_radius.get_stations_within_box(map_coords);
	stations=gps_input_pipeline.remove_blacklist(stations);
	stations.append('CME6');
	outfile_dir='Outputs/'+str(eqtime);
	return [stations, fit_table, grace_dir, start_time_infl, end_time_infl, start_time_velo, end_time_velo, N, Wn, seasonal_type, map_coords, outfile_dir];


# ------------ INPUTS  ------------- # 
def inputs(stations):
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[];
	for station_name in stations:
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, 'pbo');
		dataobj_list.append(myData);
		offsetobj_list.append(offset_obj);
		eqobj_list.append(eq_obj);
	return [dataobj_list, offsetobj_list, eqobj_list];


# ------------ COMPUTE ------------- # 
def compute(dataobj_list, offsetobj_list, eqobj_list, fit_table, grace_dir,start_time_infl, end_time_infl,start_time_velo, end_time_velo, seasonal_type, N, Wn):
	
	# Initialize output objects
	noeq_objects = []; 
	east_filt=[]; north_filt=[]; vert_filt=[]; east_inf_time=[]; north_inf_time=[]; vert_inf_time=[];
	east_change=[]; north_change=[]; vert_change=[];

	for i in range(len(dataobj_list)):
		# Remove the earthquakes and offsets
		newobj=offsets.remove_antenna_offsets(dataobj_list[i], offsetobj_list[i]);
		newobj=offsets.remove_earthquakes(newobj,eqobj_list[i]);
		newobj=gps_ts_functions.remove_outliers(newobj,15);  # 15mm horizontal outliers
		newobj=gps_seasonal_removals.make_detrended_ts(newobj, 1, seasonal_type, fit_table, grace_dir);  # can remove seasonals a few ways
		noeq_objects.append(newobj);

		print(dataobj_list[i].name);

		# Get the inflection points in the timeseries
		float_dtarray = gps_ts_functions.get_float_times(newobj.dtarray);
		[east_filtered, e_inflection_time, echange]=inflection_with_butterworth(newobj.dtarray, float_dtarray, newobj.dE, N, Wn, start_time_infl, end_time_infl,start_time_velo, end_time_velo);
		[north_filtered, n_inflection_time, nchange]=inflection_with_butterworth(newobj.dtarray, float_dtarray, newobj.dN, N, Wn, start_time_infl, end_time_infl,start_time_velo, end_time_velo);
		[vert_filtered, v_inflection_time, vchange]=inflection_with_butterworth(newobj.dtarray, float_dtarray, newobj.dU, N, Wn, start_time_infl, end_time_infl,start_time_velo, end_time_velo);

		east_filt.append(east_filtered);
		north_filt.append(north_filtered);
		vert_filt.append(vert_filtered);
		east_inf_time.append(e_inflection_time);
		north_inf_time.append(n_inflection_time);
		vert_inf_time.append(v_inflection_time);
		east_change.append(echange);
		north_change.append(nchange);
		vert_change.append(vchange);

	return [noeq_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, east_change, north_change, vert_change];


# The butterworth filter. 
def inflection_with_butterworth(dtarray, x, y, N, Wn, start_time_infl, end_time_infl,start_time_velo, end_time_velo):
	# dtarray is a datetime object array
	# x is a float
	# start_time and end_time are datetime objects
	# We search withinin starttime_infl to endtime_infl for the inflection point.
	# We use the larger array starttime_velo to endtime_velo to solve for velocity change. 

	# Build the butterworth filter
	[b,a]=butter(N, Wn, btype='low',output='ba');
	y_filtered = filtfilt(b, a, y);

	# Making a shortened filtered time series (for only the region of the earthquake)
	new_filtered_velo=[]; new_float_time_velo=[];
	new_filtered_infl=[]; new_float_time_infl=[];
	for i in range(len(dtarray)):
		if dtarray[i]>=start_time_velo and dtarray[i]<=end_time_velo:
			new_float_time_velo.append(x[i]);
			new_filtered_velo.append(y_filtered[i]);	
		if dtarray[i]>=start_time_infl and dtarray[i]<=end_time_infl:
			new_float_time_infl.append(x[i]);
			new_filtered_infl.append(y_filtered[i]);

	# Pathological case that the station has no data in the time window of interest. 
	if len(new_float_time_infl)<2:
		print("Error! Station does not have much data in the time window. Skipping.");
		return [y_filtered, dtarray[0], 0];

	# Find the inflection points. yfirst and xfirst are derivatives. 
	dy=np.diff(new_filtered_infl[:],1)
	dx=np.diff(new_float_time_infl[:],1)
	yfirst=np.divide(dy,dx);
	xfirst=np.add(new_float_time_infl[:-1],new_float_time_infl[1:])*0.5;
	slope=abs(yfirst);
	
	# Return a datetime object where the slope is minimized
	turning_point=np.argmin(slope,0);
	turning_dt=gps_ts_functions.float_to_dt(xfirst[turning_point]);

	# Get the slope of a time series
	# Find the index of the datetime where the slope went to minimum.
	turning_index_velo=new_float_time_velo.index(new_float_time_infl[turning_point]);

	if turning_index_velo<3:  # In a few rare cases, there's no good data on one side of the earthquake. 
		print("Error! Station does not have much data before the earthquake. Skipping.");
		return [y_filtered, dtarray[0], 0];
	if turning_index_velo>len(new_float_time_velo)-3:  # In a few rare cases, there's no good data on one side of the earthquake. 
		print("Error! Station does not have much after before the earthquake. Skipping.");
		return [y_filtered, dtarray[0], 0];

	slope_pre=get_ts_slope(new_float_time_velo, new_filtered_velo, 0, turning_index_velo);  # the slope before
	slope_post=get_ts_slope(new_float_time_velo, new_filtered_velo, turning_index_velo, -1);  # the slope after 
	slope_change = slope_post - slope_pre;

	return [y_filtered, turning_dt, slope_change];

def get_ts_slope(float_time, ts_data, start_index, end_index):
	slope_coef=np.polyfit(float_time[start_index:end_index],ts_data[start_index:end_index],1);
	slope=slope_coef[0];
	return slope;




# ------------ OUTPUTS ------------- # 
def outputs(noeq_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, east_change, north_change, vert_change, start_time_infl, end_time_infl, outfile_dir):
	ofile=open(outfile_dir+'_inflections.txt','w');
	for i in range(len(noeq_objects)):
		ofile.write("%s %f %f %s %s %s %.3f %.3f %.3f \n" % (noeq_objects[i].name, noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_inf_time[i], north_inf_time[i], vert_inf_time[i], east_change[i], north_change[i], vert_change[i]) );
	ofile.close();
	for i in range(len(noeq_objects)):
		output_plots(noeq_objects[i], east_filt[i], north_filt[i], vert_filt[i], east_inf_time[i], north_inf_time[i], vert_inf_time[i], start_time_infl, end_time_infl, outfile_dir);
	return;

def output_plots(noeq_obj, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, start_time, end_time, outfile_dir):

	# plt.figure(figsize=(8,8),dpi=80);
	[f,axarr]=plt.subplots(3,1,sharex=True,figsize=(8,8), dpi=80);
	axarr[0].plot_date(noeq_obj.dtarray,noeq_obj.dE,'b.');
	axarr[0].plot_date(noeq_obj.dtarray,east_filt,'r',linewidth=3);
	axarr[0].plot_date([east_inf_time, east_inf_time],[-5,5],'--k');
	axarr[0].plot_date([start_time, start_time],[-5,5],'k');
	axarr[0].plot_date([end_time, end_time],[-5,5],'k');
	axarr[0].set_ylabel('East (mm)');
	axarr[0].set_title(noeq_obj.name);
	
	axarr[1].plot_date(noeq_obj.dtarray,noeq_obj.dN,'b.');
	axarr[1].plot_date(noeq_obj.dtarray,north_filt,'r',linewidth=3);
	axarr[1].plot_date([north_inf_time, north_inf_time],[-5,5],'--k');
	axarr[1].plot_date([start_time, start_time],[-5,5],'k');
	axarr[1].plot_date([end_time, end_time],[-5,5],'k');
	axarr[1].set_ylabel('North (mm)');
	
	axarr[2].plot_date(noeq_obj.dtarray,noeq_obj.dU,'b.');
	axarr[2].plot_date(noeq_obj.dtarray,vert_filt,'r',linewidth=3);
	axarr[2].plot_date([vert_inf_time, vert_inf_time],[-15,15],'--k');
	axarr[2].plot_date([start_time, start_time],[-15,15],'k');
	axarr[2].plot_date([end_time, end_time],[-15,15],'k');
	axarr[2].set_ylabel('Up (mm)');
	axarr[2].set_xlabel('Time');
	
	plt.savefig(outfile_dir+'_'+noeq_obj.name+'_plot.png');

	return;



# --------- DRIVER ---------- # 
if __name__=="__main__":
	
	# eqtime="20140310"; starttime="20100117"; endtime="20161207"; # 2014 M6.8 Earthquake
	eqtime="20161208"; starttime="20140317"; endtime="20180901"; # 2016 M6.6 Earthquake
	# eqtime="20100110"; # # 2010 M6.5 Earthquake
	driver(eqtime, starttime, endtime);

