# Simple experiment here. How does the formal uncertainty on slope estimate vary as a function of length of the dataset, 
# or number of points in the dataset? 
# I've heard the answer is 1/n^2. Is that true? 
# Ans: no, it's n^(-3/2). 

import numpy as np 
import matplotlib.pyplot as plt 
import random
import datetime as dt
import lssq_model_errors
import gps_input_pipeline
import offsets
import gps_ts_functions
import gps_seasonal_removals


def read_real_gps_data():
	starttime=dt.datetime.strptime("20140317","%Y%m%d");
	endtime=dt.datetime.strptime("20161207","%Y%m%d");
	[newData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data('P157', 'pbo', 'NA');
	newData=offsets.remove_offsets(newData, offset_obj);
	newData=gps_ts_functions.remove_outliers(newData, 10);  # 10 mm outlier
	newData=offsets.remove_offsets(newData, eq_obj);
	newData=gps_seasonal_removals.make_detrended_ts(newData, 1, 'lssq');
	newData=gps_ts_functions.impose_time_limits(newData, starttime, endtime);
	x = gps_ts_functions.get_float_times(newData.dtarray);
	x = [i-x[0] for i in x];
	y = newData.dE;
	sig= newData.Se;
	return [x, y, sig];



[x, y, sig] = read_real_gps_data();
sig=[1*i for i in sig];
start_index=3;
num_points=[]; formal_slope_uncertainty=[]; curvefit_slope_uncertainty=[]; avr_uncertainty=[];
for i in range(start_index,len(x),10):
	avg_sigma=np.mean(sig[0:i])
	params1, covm1 = lssq_model_errors.linear_fitting_menke(x[0:i], y[0:i], avg_sigma);
	params2, covm2 = lssq_model_errors.fit_curvefit(x[0:i], y[0:i], sig[0:i]);
	params3, covm3 = lssq_model_errors.AVR(x[0:i], y[0:i], sig[0:i]);
	m_sig1=np.sqrt(covm1[0][0]);
	m_sig2=np.sqrt(covm2[0][0]);
	m_sig3=np.sqrt(covm3[0][0]);
	num_points.append(i);
	formal_slope_uncertainty.append(m_sig1);
	curvefit_slope_uncertainty.append(m_sig2);
	avr_uncertainty.append(m_sig3);



# Plotting
plt.figure();
plt.plot(num_points, formal_slope_uncertainty, '.--b',label='Formal uncertainty');
plt.plot(num_points, curvefit_slope_uncertainty, '.--r',label='Curvefit uncertainty');
plt.plot(num_points, avr_uncertainty, '.--k',label='AVR uncertainty');
y_min_uncertainty=np.multiply(2,formal_slope_uncertainty);
y_max_uncertainty=np.multiply(11, formal_slope_uncertainty);
plt.plot(num_points, y_min_uncertainty, linestyle='-',color='lightsteelblue');
plt.plot(num_points, y_max_uncertainty, linestyle='-',color='lightsteelblue');
plt.gca().fill_between(num_points, y_min_uncertainty, y_max_uncertainty, facecolor='lightsteelblue',label='2-11x formal uncertainty');
plt.xlabel('Number of Points');
plt.ylabel('Formal Uncertainty (mm/yr)');
ytest=[1000.0/(np.power(i,1.5)) for i in num_points]; # model fit
plt.plot(num_points, ytest, '--g',label='n^(-3/2) model');
plt.gca().set_yscale('log');
plt.gca().set_xscale('log');
plt.grid(True);
plt.legend();
plt.savefig('formal_uncertainty.eps');
