# General tools used for tremor analysis.


import numpy as np 
import datetime as dt 
import collections
import gps_input_pipeline
import offsets
import gps_ts_functions
import gps_seasonal_removals

TremorCat = collections.namedtuple("TremorCat",['dtarray','lonarray','latarray']);


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

	