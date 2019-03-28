# Python script
# March 27, 2019
# The goal of this script is to test a variety of seasonal sinusoids with reasonable amplitude and random phase, 
# And see if we get the 2014-2016 velocity changes. 

import numpy as np
import matplotlib.pyplot as plt 
import datetime as dt 
import collections
import gps_ts_functions
import gps_input_pipeline
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

def configure():
	amplitude=1;  # mm, one-half the peak-to-peak amplitude. 
	EQ1=dt.datetime.strptime("20100703","%Y%m%d");  # earthquake was 20100110. We want to skip six months. 
	EQ2=dt.datetime.strptime("20140310","%Y%m%d");
	EQ3=dt.datetime.strptime("20161208","%Y%m%d");
	EQ4=dt.datetime.strptime("20180915","%Y%m%d");
	return amplitude, EQ1, EQ2, EQ3, EQ4; 

# Basic idea: 
# Amplitude^2 = A^2 + B^2
def compute(gps_data, amplitude, EQ1, EQ2, EQ3, EQ4):
	# Setting up
	decyear=gps_ts_functions.get_float_times(gps_data.dtarray);
	blanks=np.zeros(np.shape(decyear));
	time_after_start_date=7;
	start_time=dt.datetime.strptime("20100603","%Y%m%d");
	end_time=dt.datetime.strptime("20180930","%Y%m%d");
	A_coef=np.arange(-1,1,0.04);
	A_axis=[];
	dv_2014=[];
	dv_2016=[];

	# Starting to make a complicated figure with a sub-plot inside the middle. 
	fig = plt.figure(figsize=(10,10));	
	ax1 = fig.add_axes([0.1,0.1,0.8,0.8]);
	ax2 = fig.add_axes([0.52,0.7,0.35,0.1]);


	for A in A_coef:
		B=np.sqrt(amplitude*amplitude - A*A);
		myfunction = gps_ts_functions.annual_only_function(decyear,[A,B]);
		myts = Timeseries(name='TEST',coords=[0,0],dtarray=gps_data.dtarray, dN=blanks, dE=myfunction, dU=blanks, Sn=blanks, Se=blanks, Su=blanks, EQtimes=[]);
		[T2slope, _, _, _, _, _]=gps_ts_functions.get_slope(myts,starttime=EQ1+dt.timedelta(days=time_after_start_date),endtime=EQ2);
		[T3slope, _, _, _, _, _]=gps_ts_functions.get_slope(myts,starttime=EQ2+dt.timedelta(days=time_after_start_date),endtime=EQ3);
		[T4slope, _, _, _, _, _]=gps_ts_functions.get_slope(myts,starttime=EQ3+dt.timedelta(days=time_after_start_date),endtime=EQ4);
		A_axis.append(A);
		dv_2014.append(T3slope-T2slope);
		dv_2016.append(T4slope-T3slope);

		B=-B;
		myfunction = gps_ts_functions.annual_only_function(decyear,[A,B]);
		myts = Timeseries(name='TEST',coords=[0,0],dtarray=gps_data.dtarray, dN=blanks, dE=myfunction, dU=blanks, Sn=blanks, Se=blanks, Su=blanks, EQtimes=[]);
		[T2slope, _, _, _, _, _]=gps_ts_functions.get_slope(myts,starttime=EQ1+dt.timedelta(days=time_after_start_date),endtime=EQ2);
		[T3slope, _, _, _, _, _]=gps_ts_functions.get_slope(myts,starttime=EQ2+dt.timedelta(days=time_after_start_date),endtime=EQ3);
		[T4slope, _, _, _, _, _]=gps_ts_functions.get_slope(myts,starttime=EQ3+dt.timedelta(days=time_after_start_date),endtime=EQ4);
		A_axis.append(A);
		dv_2014.append(T3slope-T2slope);
		dv_2016.append(T4slope-T3slope);	

		ax2.plot(myts.dtarray, myfunction);  
		# I want to put a color mapable here for the phase. 

	ax2.plot([EQ1,EQ1],[-amplitude,amplitude],'--k');
	ax2.plot([EQ2,EQ2],[-amplitude,amplitude],'--k');
	ax2.plot([EQ3,EQ3],[-amplitude,amplitude],'--k');
	ax2.plot([EQ4,EQ4],[-amplitude,amplitude],'--k');
	ax2.text(dt.datetime.strptime("20120106","%Y%m%d"),1.2,"T2",fontsize=16,fontweight='bold');
	ax2.text(dt.datetime.strptime("20150606","%Y%m%d"),1.2,"T3",fontsize=16,fontweight='bold');
	ax2.text(dt.datetime.strptime("20171006","%Y%m%d"),1.2,"T4",fontsize=16,fontweight='bold');
	ax2.set_xlim(start_time, end_time);

	ax1.plot(A_axis, dv_2014, '.', label='Apparent East Delta-V, T3-T2');
	ax1.plot(A_axis, dv_2016, '.', label='Apparent East Delta-V, T4-T3');
	ax1.plot([-1,1],[amplitude, amplitude],'--k');
	ax1.text(-1,amplitude-0.08,'Imposed Seasonal Amplitude',fontsize=16);
	ax1.set_ylabel('Apparent DeltaV (mm/yr)',fontsize=16);
	ax1.set_xlabel('Random Phase',fontsize=16);
	ax1.grid('on');
	ax1.legend(fontsize=16);
	plt.savefig('Synethetic_Phase_deltaV.eps');
	return A_axis, dv_2014, dv_2016; 


if __name__=="__main__":
	amplitude, EQ1, EQ2, EQ3, EQ4 = configure();
	[gps_ts,_,_] = gps_input_pipeline.get_station_data("P160",'pbo');
	A_coef, dv_2014, dv_2016 = compute(gps_ts, amplitude, EQ1, EQ2, EQ3, EQ4);



