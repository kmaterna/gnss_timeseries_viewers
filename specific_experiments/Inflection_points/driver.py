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


# Reference : 
# Timeseries = collections.namedtuple("Timeseries",
#      ['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);


def driver(eqtime):
	[stations, filenames, earthquakes_dir, offsets_dir, start_time, end_time, N, Wn, avoid_edges, map_coords, outfile_dir] = configure(eqtime);
	dataobj_list = inputs(stations, filenames);
	[noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time]=compute(dataobj_list, earthquakes_dir, offsets_dir, start_time, end_time, N, Wn, avoid_edges);
	outputs(noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, outfile_dir);
	return;


# ------------ CONFIGURE ------------- # 
def configure(eqtime):
	EQtime  = dt.datetime.strptime(eqtime, "%Y%m%d");
	pbo_velfile="../../GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt";
	earthquakes_dir = "../../GPS_POS_DATA/PBO_Event_Files/";
	offsets_dir = "../../GPS_POS_DATA/Offsets/";
	pre_event_duration = 2.5; # years
	post_event_duration = 2.5; # years
	start_time=EQtime-dt.timedelta(days=pre_event_duration*365);
	end_time=EQtime+dt.timedelta(days=post_event_duration*365);

	# Butterworth parameters
	N=3;  # Order of butterworth filter
	Wn=1/365.0;  # 1/period (days) of cutoff frequency. 	
	avoid_edges=730;  # When you find inflection points, you avoid each end of the filter (days)
	print('searching %f years, centered on the event time' % ( ( (pre_event_duration+post_event_duration)*365.25-2*avoid_edges)/365.24 ) );
	# Presently, I'm searching 1 year surrounding the event time for an inflection point. 

	map_coords=[-125, -122, 39, 41.5];
	stations = stations_within_radius.get_stations_within_box(map_coords,pbo_velfile);
	# stations=['P659'];
	filenames=[];
	for station in stations:
		filenames.append("../../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos");
	outfile_dir='Outputs/'+str(eqtime);
	return [stations, filenames, earthquakes_dir, offsets_dir, start_time, end_time, N, Wn, avoid_edges, map_coords, outfile_dir];


# ------------ INPUTS  ------------- # 
def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;


# ------------ COMPUTE ------------- # 
def compute(dataobj_list, earthquakes_dir, offsets_dir, start_time, end_time, N, Wn, avoid_edges):
	
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
		[east_filtered, e_inflection_time]=inflection_with_butterworth(float_dtarray, short_obj.dE, N, Wn, avoid_edges);
		[north_filtered, n_inflection_time]=inflection_with_butterworth(float_dtarray, short_obj.dN, N, Wn, avoid_edges);
		[vert_filtered, v_inflection_time]=inflection_with_butterworth(float_dtarray, short_obj.dU, N, Wn, avoid_edges);

		east_filt.append(east_filtered);
		north_filt.append(north_filtered);
		vert_filt.append(vert_filtered);
		east_inf_time.append(e_inflection_time);
		north_inf_time.append(n_inflection_time);
		vert_inf_time.append(v_inflection_time);

	return [noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time];


# The butterworth filter. 
def inflection_with_butterworth(x, y, N, Wn, avoid_edges):

	# Build the butterworth filter
	[b,a]=butter(N, Wn, btype='low',output='ba');
	y_filtered = filtfilt(b, a, y);

	# Avoid the edges of the filter. 
	if len(y_filtered)<2*avoid_edges:
		print("Error: avoid_edges is too long. Cutting to a smaller value. ");
		avoid_edges=int((len(y_filtered)-365)/2);

	# Find the inflection points. yfirst and xfirst are derivatives. 
	dy=np.diff(y_filtered[avoid_edges:-avoid_edges],1)
	dx=np.diff(x[avoid_edges:-avoid_edges],1)
	yfirst=np.divide(dy,dx);
	xfirst=np.add(x[avoid_edges:-avoid_edges-1],x[avoid_edges+1:-avoid_edges])*0.5;

	# Return a datetime object where the slope is minimized
	slope=abs(yfirst);
	turning_point=np.argmin(slope,0);
	turning_dt=gps_ts_functions.float_to_dt(xfirst[turning_point]); 

	return [y_filtered, turning_dt];





# ------------ OUTPUTS ------------- # 
def outputs(noeq_objects, short_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, outfile_dir):
	ofile=open(outfile_dir+'_inflections.txt','w');
	for i in range(len(noeq_objects)):
		ofile.write("%f %f %s %s %s\n" % (noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_inf_time[i], north_inf_time[i], vert_inf_time[i]) );
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
	axarr[0].plot_date([east_inf_time, east_inf_time],[-5,5],'--k');
	axarr[0].set_ylabel('East (mm)');
	axarr[0].set_title(noeq_obj.name);
	
	axarr[1].plot_date(noeq_obj.dtarray,noeq_obj.dN);
	axarr[1].plot_date(short_obj.dtarray,short_obj.dN,'r.');
	axarr[1].plot_date(short_obj.dtarray,north_filt,'g',linewidth=3);
	axarr[1].plot_date([north_inf_time, north_inf_time],[-5,5],'--k');
	axarr[1].set_ylabel('North (mm)');
	
	axarr[2].plot_date(noeq_obj.dtarray,noeq_obj.dU);
	axarr[2].plot_date(short_obj.dtarray,short_obj.dU,'r.');
	axarr[2].plot_date(short_obj.dtarray,vert_filt,'g',linewidth=3);
	axarr[2].plot_date([vert_inf_time, vert_inf_time],[-15,15],'--k');
	axarr[2].set_ylabel('Up (mm)');
	axarr[2].set_xlabel('Time');
	
	plt.savefig(outfile_dir+'_'+noeq_obj.name+'_plot.png');

	return;



# --------- DRIVER ---------- # 
if __name__=="__main__":
	
	# eqtime="20140314"; # 2014 M6.8 Earthquake
	eqtime="20161208"; # # 2016 M6.6 Earthquake
	# eqtime="20100110"; # # 2010 M6.5 Earthquake
	driver(eqtime);

