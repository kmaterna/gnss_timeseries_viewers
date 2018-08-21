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
	[stations, filenames, earthquakes_dir, offsets_dir, fit_table, grace_dir,start_time, end_time, N, Wn, seasonal_type, map_coords, outfile_dir] = configure(eqtime);
	dataobj_list = inputs(stations, filenames);
	[noeq_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, east_change, north_change, vert_change]=compute(dataobj_list, earthquakes_dir, offsets_dir, fit_table, grace_dir,start_time, end_time, seasonal_type, N, Wn);
	outputs(noeq_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, east_change, north_change, vert_change, start_time, end_time, outfile_dir);
	return;


# ------------ CONFIGURE ------------- # 
def configure(eqtime):
	EQtime  = dt.datetime.strptime(eqtime, "%Y%m%d");
	pbo_velfile="../../GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt";
	earthquakes_dir = "../../GPS_POS_DATA/PBO_Event_Files/";
	offsets_dir = "../../GPS_POS_DATA/Offsets/";
	fit_table="../../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt"
	grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"
	seasonal_type="notch";

	if eqtime[0:4]=='2016' or eqtime[0:4]=='2017':
		pre_event_duration = 0.5; # years
		post_event_duration = 0.5; # years
		print('searching %f years, centered on the event time' % ( pre_event_duration+post_event_duration ) );
		# Presently, I'm searching 1 year surrounding the event time for an inflection point. 	
	else:
		pre_event_duration = 0.5; # years
		post_event_duration = 0.5; # years
		print('searching %f years, centered on the event time' % ( pre_event_duration+post_event_duration) );

	start_time=EQtime-dt.timedelta(days=pre_event_duration*365);
	end_time=EQtime+dt.timedelta(days=post_event_duration*365);

	# Butterworth parameters
	N=3;  # Order of butterworth filter
	Wn=1/365.0;  # 1/period (days) of cutoff frequency. 	

	map_coords=[-125, -122, 39, 41.5];
	stations = stations_within_radius.get_stations_within_box(map_coords,pbo_velfile);
	# stations=['P659'];
	filenames=[];
	for station in stations:
		filenames.append("../../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos");
	outfile_dir='Outputs/'+str(eqtime);
	return [stations, filenames, earthquakes_dir, offsets_dir,fit_table, grace_dir, start_time, end_time, N, Wn, seasonal_type, map_coords, outfile_dir];


# ------------ INPUTS  ------------- # 
def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;


# ------------ COMPUTE ------------- # 
def compute(dataobj_list, earthquakes_dir, offsets_dir, fit_table, grace_dir,start_time, end_time, seasonal_type, N, Wn):
	
	# Initialize output objects
	noeq_objects = []; 
	east_filt=[]; north_filt=[]; vert_filt=[]; east_inf_time=[]; north_inf_time=[]; vert_inf_time=[];
	east_change=[]; north_change=[]; vert_change=[];

	for i in range(len(dataobj_list)):
		# Remove the earthquakes
		newobj=gps_ts_functions.remove_offsets(dataobj_list[i],offsets_dir);
		newobj=gps_ts_functions.remove_earthquakes(newobj,earthquakes_dir);
		newobj=gps_ts_functions.remove_outliers(newobj,15);  # 15mm horizontal outliers
		newobj=gps_ts_functions.make_detrended_option(newobj, 1, seasonal_type, fit_table, grace_dir);  # can remove seasonals a few ways
		noeq_objects.append(newobj);

		# Get the inflection points in the timeseries
		float_dtarray = gps_ts_functions.get_float_times(newobj.dtarray);
		[east_filtered, e_inflection_time, echange]=inflection_with_butterworth(newobj.dtarray, float_dtarray, newobj.dE, N, Wn, start_time, end_time);
		[north_filtered, n_inflection_time, nchange]=inflection_with_butterworth(newobj.dtarray, float_dtarray, newobj.dN, N, Wn, start_time, end_time);
		[vert_filtered, v_inflection_time, vchange]=inflection_with_butterworth(newobj.dtarray, float_dtarray, newobj.dU, N, Wn, start_time, end_time);

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
def inflection_with_butterworth(dtarray, x, y, N, Wn, start_time, end_time):

	# Build the butterworth filter
	[b,a]=butter(N, Wn, btype='low',output='ba');
	y_filtered = filtfilt(b, a, y);

	# Making a shortened filtered time series (for only the region of the earthquake)
	new_filtered=[]; new_float_time=[];
	for i in range(len(dtarray)):
		if dtarray[i]>=start_time and dtarray[i]<=end_time:
			new_float_time.append(x[i]);
			new_filtered.append(y_filtered[i]);			

	# Find the inflection points. yfirst and xfirst are derivatives. 
	dy=np.diff(new_filtered[:],1)
	dx=np.diff(new_float_time[:],1)
	yfirst=np.divide(dy,dx);
	xfirst=np.add(new_float_time[:-1],new_float_time[1:])*0.5;

	# Return a datetime object where the slope is minimized
	slope=abs(yfirst);
	turning_point=np.argmin(slope,0);
	turning_dt=gps_ts_functions.float_to_dt(xfirst[turning_point]);

	# Get the slope of a time series, but make sure that you can use at least a month on either side. 
	if turning_point<30 or turning_point>len(yfirst)-30:
		slope_change=0;
	else:
		slope_pre=get_ts_slope(new_float_time, new_filtered, 0, turning_point);  # the slope before
		slope_post=get_ts_slope(new_float_time, new_filtered, turning_point, -1);  # the slope after 
		slope_change = slope_post - slope_pre;

	return [y_filtered, turning_dt, slope_change];

def get_ts_slope(float_time, ts_data, start_index, end_index):
	slope_coef=np.polyfit(float_time[start_index:end_index],ts_data[start_index:end_index],1);
	slope=slope_coef[0];
	return slope;




# ------------ OUTPUTS ------------- # 
def outputs(noeq_objects, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, east_change, north_change, vert_change, start_time, end_time, outfile_dir):
	ofile=open(outfile_dir+'_inflections.txt','w');
	for i in range(len(noeq_objects)):
		ofile.write("%s %f %f %s %s %s %.3f %.3f %.3f \n" % (noeq_objects[i].name, noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_inf_time[i], north_inf_time[i], vert_inf_time[i], east_change[i], north_change[i], vert_change[i]) );
	ofile.close();
	for i in range(len(noeq_objects)):
		output_plots(noeq_objects[i], east_filt[i], north_filt[i], vert_filt[i], east_inf_time[i], north_inf_time[i], vert_inf_time[i], start_time, end_time, outfile_dir);
	return;

def output_plots(noeq_obj, east_filt, north_filt, vert_filt, east_inf_time, north_inf_time, vert_inf_time, start_time, end_time, outfile_dir):

	plt.figure();
	[f,axarr]=plt.subplots(3,1,sharex=True);
	axarr[0].plot_date(noeq_obj.dtarray,noeq_obj.dE);
	axarr[0].plot_date(noeq_obj.dtarray,east_filt,'g',linewidth=3);
	axarr[0].plot_date([east_inf_time, east_inf_time],[-5,5],'--k');
	axarr[0].plot_date([start_time, start_time],[-5,5],'k');
	axarr[0].plot_date([end_time, end_time],[-5,5],'k');
	axarr[0].set_ylabel('East (mm)');
	axarr[0].set_title(noeq_obj.name);
	
	axarr[1].plot_date(noeq_obj.dtarray,noeq_obj.dN);
	axarr[1].plot_date(noeq_obj.dtarray,north_filt,'g',linewidth=3);
	axarr[1].plot_date([north_inf_time, north_inf_time],[-5,5],'--k');
	axarr[1].plot_date([start_time, start_time],[-5,5],'k');
	axarr[1].plot_date([end_time, end_time],[-5,5],'k');
	axarr[1].set_ylabel('North (mm)');
	
	axarr[2].plot_date(noeq_obj.dtarray,noeq_obj.dU);
	axarr[2].plot_date(noeq_obj.dtarray,vert_filt,'g',linewidth=3);
	axarr[2].plot_date([vert_inf_time, vert_inf_time],[-15,15],'--k');
	axarr[2].plot_date([start_time, start_time],[-15,15],'k');
	axarr[2].plot_date([end_time, end_time],[-15,15],'k');
	axarr[2].set_ylabel('Up (mm)');
	axarr[2].set_xlabel('Time');
	
	plt.savefig(outfile_dir+'_'+noeq_obj.name+'_plot.png');

	return;



# --------- DRIVER ---------- # 
if __name__=="__main__":
	
	eqtime="20140314"; # 2014 M6.8 Earthquake
	# eqtime="20161208"; # # 2016 M6.6 Earthquake
	# eqtime="20100110"; # # 2010 M6.5 Earthquake
	driver(eqtime);

