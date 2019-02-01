# This is the driver for a new program 

import numpy as np 
import matplotlib.pyplot as plt 
import random
import datetime as dt
import lssq_model_errors
import gps_input_pipeline
import offsets
import gps_ts_functions
import gps_seasonal_removals



def generate_data_const_sig():
	# Test data, such as deformation vs. time 
	x=np.linspace(0,100,10);

	# True parameters:
	m=1;
	b=0;
	dsig=20;

	# The linear model
	y=[m*i + b + dsig*(-0.5+random.random()) for i in x];
	sig=[dsig for i in x];
	return [x, y, sig];

def read_real_gps_data():
	starttime=dt.datetime.strptime("20140317","%Y%m%d");
	endtime=dt.datetime.strptime("20161207","%Y%m%d");
	# starttime=dt.datetime.strptime("20161210","%Y%m%d");
	# endtime=dt.datetime.strptime("20180915","%Y%m%d");	
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

def plot_data(x,y,sig,params1, covm1, params2, covm2, params3, covm3):
	f=plt.figure();
	plt.errorbar(x,y,yerr=sig,fmt='o',zorder=0);
	
	# Plotting option 1
	m_model=params1[0]; b_model=params1[1]; 
	m_sig=np.sqrt(covm1[0][0]);
	y_model=produce_linear_model(x, m_model, b_model);
	y_max=produce_linear_model(x, m_model+m_sig, b_model);
	y_min=produce_linear_model(x, m_model-m_sig, b_model);
	plt.plot(x,y_model,'k');
	plt.plot(x,y_max,'--k');
	p1 = plt.plot(x,y_min,'--k',label='Formal Unc. Prop.');

	# Plotting option 2
	m_model=params2[0]; b_model=params2[1]; 
	m_sig=np.sqrt(covm2[0][0]);
	y_model=produce_linear_model(x, m_model, b_model);
	y_max=produce_linear_model(x, m_model+m_sig, b_model);
	y_min=produce_linear_model(x, m_model-m_sig, b_model);
	plt.plot(x,y_model,'r');
	plt.plot(x,y_max,'--r');
	p2=plt.plot(x,y_min,'--r',label='scipy curvefit');

	# Plotting option 3
	m_model=params3[0]; b_model=params3[1]; 
	m_sig=np.sqrt(covm3[0][0]);
	y_model=produce_linear_model(x, m_model, b_model);
	y_max=produce_linear_model(x, m_model+m_sig, b_model);
	y_min=produce_linear_model(x, m_model-m_sig, b_model);
	plt.plot(x,y_model,color='purple');
	plt.plot(x,y_max,linestyle='--',color='purple');
	p2=plt.plot(x,y_min,linestyle='--',color='purple',label='AVR');
	plt.xlabel('Time (years)');
	plt.ylabel('Position (mm)');

	plt.legend();

	plt.savefig('model.eps');
	return;

def produce_linear_model(x, m, b):
	y=[i*m+b for i in x];
	return y;


if __name__=="__main__":
	# [x, y, sig] = generate_data_const_sig();
	[x, y, sig] = read_real_gps_data();
	avg_sigma = np.mean(sig);
	params1, covm1 = lssq_model_errors.linear_fitting_menke(x, y, avg_sigma);
	params2, covm2 = lssq_model_errors.fit_curvefit(x, y, sig);
	params3, covm3 = lssq_model_errors.AVR(x, y, sig);
	plot_data(x,y,sig, params1, covm1, params2, covm2, params3, covm3);

