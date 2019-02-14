# General tools used for tremor analysis.


import numpy as np 
import datetime as dt 
import collections
import gps_input_pipeline
import offsets
import gps_ts_functions
import gps_seasonal_removals

TremorCat = collections.namedtuple("TremorCat",['dtarray','lonarray','latarray']);
TremorCatDepths = collections.namedtuple("TremorCatDepths",['dtarray','lonarray','latarray','depth']);


def restrict_to_box(tremor, box_interest, start_time, end_time):
	newdt=[]; newlon=[]; newlat=[];
	for i in range(len(tremor.dtarray)):
		if tremor.dtarray[i]>start_time and tremor.dtarray[i]<end_time:
			if tremor.lonarray[i]>box_interest[0] and tremor.lonarray[i]<box_interest[1]:
				if tremor.latarray[i]>box_interest[2] and tremor.latarray[i]<box_interest[3]:
					newdt.append(tremor.dtarray[i]);
					newlon.append(tremor.lonarray[i]);
					newlat.append(tremor.latarray[i]);
	newtremor=TremorCat(dtarray=newdt, lonarray=newlon, latarray=newlat);
	return newtremor;


def get_cumulative_plot(tremor, box_interest, start_time, end_time):
	# Returns two arrays that can be plotted against each other to give the cumulative tremor plot. 
	dt_interest=[];
	cnumber=[];
	dt_interest.append(start_time);
	cnumber.append(0);
	for i in range(len(tremor.dtarray)):
		if tremor.dtarray[i]>start_time:
			if tremor.lonarray[i]>box_interest[0] and tremor.lonarray[i]<box_interest[1]:
				if tremor.latarray[i]>box_interest[2] and tremor.latarray[i]<box_interest[3]:
					dt_interest.append(tremor.dtarray[i]);
					dt_interest.append(tremor.dtarray[i]);
					cnumber.append(cnumber[-1]);
					cnumber.append(cnumber[-1]+1);
	cnumber=np.array(cnumber);	
	return [dt_interest, cnumber];



def get_rates(tremor):
	# Returns two arrays that can be potted against each other to give the tremor rate. 
	# Don't give this too much data. 
	dts=[];
	rate=[];
	dts.append(tremor.dtarray[0]);
	rate.append(0);

	while dts[-1]<tremor.dtarray[-1]:
		dts.append(dts[-1]+dt.timedelta(days=1));

		counts_that_day=0;
		for i in range(len(tremor.dtarray)):
			if tremor.dtarray[i]>dts[-1]:
				break;
			if tremor.dtarray[i]<dts[-2]:
				continue;
			if tremor.dtarray[i]>dts[-2] and tremor.dtarray[i]<dts[-1]:
				counts_that_day=counts_that_day+1;
		rate.append(counts_that_day);

	return [dts, rate]; 


def get_detrended_gps_station(station_name):
	datasource='pbo';
	[newData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, datasource);
	newData=offsets.remove_antenna_offsets(newData, offset_obj);
	newData=gps_ts_functions.remove_outliers(newData, 5);  # mm for outliers
	newData=offsets.remove_earthquakes(newData, eq_obj);
	trend_out=gps_seasonal_removals.make_detrended_ts(newData, 1, 'notch');
	return trend_out;


def write_deltat_certain_day(tremor, eqtime, endtime, outfile,unit='hours'):
	# This is for GMT plotting. If we want to write the tremor within a time window of a certain event.
	# Will write lon, lat, time delay in hours
	ofile=open(outfile,'w');
	for i in range(len(tremor.dtarray)):
		if tremor.dtarray[i]>=eqtime:
			if tremor.dtarray[i]<=endtime:
				deltat=tremor.dtarray[i]-eqtime;
				if unit=='hours':
					ofile.write("%f %f %f \n" % (tremor.lonarray[i], tremor.latarray[i],deltat.days*24+deltat.seconds/3600.0) ); # in hours
				else:
					ofile.write("%f %f %f \n" % (tremor.lonarray[i], tremor.latarray[i],deltat.days+deltat.seconds/(24*3600.0)) );  # in days
	ofile.close();
	return;
	


def get_depth_projection(coords, xdata, ydata, zdata):
	# Coords : a list of two-valued arrays
	# xdata  : an ascending array of x-values from the fault grd file
	# ydata  : an ascending array of y-values from the fault grd file
	# xdata  : a 2D array of depths from the fault grd file
	# Returns the depths of the tremor locations on the underlying xyz fault surface
	result=[];
	for i in range(len(coords)):
		# Defensive programming
		if coords[i][0]>np.max(xdata) or coords[i][0]<np.min(xdata):
			result.append(np.nan);
			continue;
		if coords[i][1]>np.max(ydata) or coords[i][1]<np.min(ydata):
			result.append(np.nan);
			continue;
		
		# Find the part of the array where the given data is located. 
		idx_left=np.searchsorted(xdata, coords[i][0], side="left");  # the lower value
		idx_right=np.searchsorted(xdata, coords[i][0], side="right");  # the higher value
		idy_left=np.searchsorted(ydata, coords[i][1], side="left");  # the lower value
		idy_right=np.searchsorted(ydata, coords[i][1], side="right");  # the higher value
		x  = coords[i][0];   # the x coordinate we're going to interpolate over. 
		xp = [idx_left, idx_right];
		fp_left  = np.mean([zdata[idy_left][idx_left], zdata[idy_right][idx_left]]);
		fp_right = np.mean([zdata[idy_left][idx_right], zdata[idy_right][idx_right]]);
		fp = [fp_left, fp_right];
		depth = np.interp(x, xp, fp);
		result.append(depth);
	return result;


def associate_depths(tremor, depths):
	tremor_with_depths = TremorCatDepths(dtarray=tremor.dtarray, lonarray=tremor.lonarray, latarray=tremor.latarray, depth=depths);
	return tremor_with_depths;





