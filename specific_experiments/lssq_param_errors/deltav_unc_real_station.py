# A few lines of code to get a simple velocity uncertainty estimate for a real station 
# Prototype for what's going into the larger velocity difference estimator. 

import numpy as np 
import datetime as dt
import lssq_model_errors
import gps_input_pipeline
import offsets
import gps_ts_functions
import gps_seasonal_removals

def read_real_gps_data(station, starttime, endtime):
	[newData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station, 'pbo',, 'NA';
	newData=offsets.remove_offsets(newData, offset_obj);
	newData=gps_ts_functions.remove_outliers(newData, 10);  # 10 mm outlier
	newData=offsets.remove_offsets(newData, eq_obj);
	newData=gps_seasonal_removals.make_detrended_ts(newData, 1, 'lssq');
	newData=gps_ts_functions.impose_time_limits(newData, starttime, endtime);
	x = gps_ts_functions.get_float_times(newData.dtarray);
	y = newData.dE;
	sig= newData.Se;
	return [x, y, sig];

def get_slope_uncertainty(station, starttime, endtime):
	[x, y, sig] = read_real_gps_data(station, starttime, endtime);
	params, covm = lssq_model_errors.AVR(x, y, sig);
	slope = params[0];
	sigma = np.sqrt(covm[0][0]);
	return slope, sigma;

if __name__=="__main__":
	station = "P157";
	starttime1=dt.datetime.strptime("20140317","%Y%m%d");
	endtime1=dt.datetime.strptime("20161207","%Y%m%d");	
	starttime2=dt.datetime.strptime("20161210","%Y%m%d");
	endtime2=dt.datetime.strptime("20181210","%Y%m%d");

	slope1, sigma1 = get_slope_uncertainty(station, starttime1, endtime1);
	slope2, sigma2 = get_slope_uncertainty(station, starttime2, endtime2);
	deltaV = abs(slope1-slope2);
	overall_sigma = sigma1+sigma2;
	print("Station "+station);
	print("Interval 1: %f +- %f mm/yr" % (slope1, sigma1) );
	print("Interval 2: %f +- %f mm/yr" % (slope2, sigma2) );
	print("Delta V : %f +- %f" % (deltaV, overall_sigma) );
