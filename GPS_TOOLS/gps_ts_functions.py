# May/June 2018
# This is a toolbox that operates on Timeseries objects. 
# Contains functions to map, filter, and reduce generic GPS time series
# 
# FIRST TYPE OF FUNCTION: RETURNS TSOBJECT
# remove_outliers(Data0, outliers_def);
# impose_time_limits(Data0, starttime, endtime)
# remove_nans(Data0)
# detrend_data_by_value(Data0,east_params,north_params,vert_params)
# remove_seasonal_by_value(Data0,east_params,north_params,vert_params)
# pair_gps_model(gps_data, model_data)
# pair_gps_model_keeping_gps(gps_data, model_data)
# get_referenced_data(roving_station_data, base_station_data)
# remove_constant(Data0, east_offset, north_offset, vert_offset)
# 
# SECOND TYPE OF FUNCTION: RETURNS VALUES
# get_slope(Data0, starttime=[], endtime=[],missing_fraction=0.6)
# get_slope_unc(dataObj, starttime, endtime)
# get_linear_annual_semiannual(Data0, starttime=[], endtime=[],critical_len=365)
# get_means(Data0, starttime=[], endtime=[])
# get_values_at_date(Data0, selected_date)

import numpy as np 
import subprocess, sys, collections
import datetime as dt 
from scipy import signal
import gps_io_functions
import lssq_model_errors 
from scipy.optimize import curve_fit

# A line for referencing the namedtuple definition. 
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


# -------------------------------------------- # 
# FUNCTIONS THAT TAKE TIME SERIES OBJECTS # 
# AND RETURN OTHER TIME SERIES OBJECTS # 
# -------------------------------------------- # 


def remove_outliers(Data0, outliers_def):
	medfilt_e=signal.medfilt(Data0.dE, 35);
	medfilt_n=signal.medfilt(Data0.dN, 35);
	medfilt_u=signal.medfilt(Data0.dU, 35);
	newdt=[]; newdN=[]; newdE=[]; newdU=[]; newSe=[]; newSn=[]; newSu=[]; 
	for i in range(len(medfilt_e)):
		if abs(Data0.dE[i]-medfilt_e[i])<outliers_def and abs(Data0.dN[i]-medfilt_n[i])<outliers_def and abs(Data0.dU[i]-medfilt_u[i])<outliers_def*2:
			newdt.append(Data0.dtarray[i]);
			newdE.append(Data0.dE[i]);
			newdN.append(Data0.dN[i]);
			newdU.append(Data0.dU[i]);
			newSe.append(Data0.Se[i]);
			newSn.append(Data0.Sn[i]);
			newSu.append(Data0.Su[i]);
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdt, dN=newdN, dE=newdE, dU=newdU, Sn=newSn, Se=newSe, Su=newSu, EQtimes=Data0.EQtimes);
	return newData;


def impose_time_limits(Data0, starttime, endtime):
	# Starttime and endtime are datetime objects
	newdtarray=[]; newdN=[]; newdE=[]; newdU=[]; newSn=[]; newSe=[]; newSu=[];
	for i in range(len(Data0.dN)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime:
			newdtarray.append(Data0.dtarray[i]);
			newdE.append(Data0.dE[i]);
			newdN.append(Data0.dN[i]);
			newdU.append(Data0.dU[i]);
			newSe.append(Data0.Se[i]);
			newSn.append(Data0.Sn[i]);
			newSu.append(Data0.Su[i]);
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=newSn, Se=newSe, Su=newSu, EQtimes=Data0.EQtimes);
	return newData;


def remove_nans(Data0):
	idxE=np.isnan(Data0.dE);
	idxN=np.isnan(Data0.dN);
	idxU=np.isnan(Data0.dU);
	temp_dates=[];
	temp_east=[];
	temp_north=[];
	temp_vert=[];
	temp_Sn=[];
	temp_Se=[];
	temp_Su=[];
	for i in range(len(Data0.dtarray)):
		if idxE[i]==0 and idxN[i]==0 and idxU[i]==0:
			temp_dates.append(Data0.dtarray[i]);
			temp_east.append(Data0.dE[i]);
			temp_north.append(Data0.dN[i]);
			temp_vert.append(Data0.dU[i]);
			temp_Se.append(Data0.Se[i]);
			temp_Sn.append(Data0.Sn[i]);
			temp_Su.append(Data0.Su[i]);
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=temp_dates, dN=temp_north, dE=temp_east, dU=temp_vert, Sn=temp_Sn, Se=temp_Se, Su=temp_Su, EQtimes=Data0.EQtimes);
	return newData;

def detrend_data_by_value(Data0,east_params,north_params,vert_params):
	if sum(np.isnan(east_params))>0 or sum(np.isnan(north_params))>0 or sum(np.isnan(vert_params))>0:
		print("ERROR: Your input slope values contain nan!");
		return Data0;

	# Parameters Format: slope, a2(cos), a1(sin), s2, s1. 
	east_detrended=[]; north_detrended=[]; vert_detrended=[];
	idx=np.isnan(Data0.dE);
	if(sum(idx))>0:  # if there are nans, please pull them out. 
		Data0=remove_nans(Data0);
	decyear=get_float_times(Data0.dtarray);
	
	east_model=linear_annual_semiannual_function(decyear,east_params);
	north_model=linear_annual_semiannual_function(decyear,north_params);
	vert_model=linear_annual_semiannual_function(decyear,vert_params);

	for i in range(len(decyear)):
		east_detrended.append(Data0.dE[i]-(east_model[i]) );
		north_detrended.append(Data0.dN[i]-(north_model[i]) );
		vert_detrended.append(Data0.dU[i]-(vert_model[i]) );
	east_detrended=east_detrended-east_detrended[0];
	north_detrended=north_detrended-north_detrended[0];
	vert_detrended=vert_detrended-vert_detrended[0];
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=north_detrended, dE=east_detrended, dU=vert_detrended, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
	return newData;

def remove_seasonal_by_value(Data0,east_params,north_params,vert_params):
	# Least squares seasonal parameters. Remove seasonal components. 
	# Parameters Format: slope, a2(cos), a1(sin), s2, s1. 
	east_detrended=[]; north_detrended=[]; vert_detrended=[];
	idx=np.isnan(Data0.dE);
	if(sum(idx))>0:  # if there are nans, please pull them out. 
		Data0=remove_nans(Data0);
	decyear=get_float_times(Data0.dtarray);
	
	east_model=	annual_semiannual_only_function(decyear,east_params[1:]);
	north_model=annual_semiannual_only_function(decyear,north_params[1:]);
	vert_model=annual_semiannual_only_function(decyear,vert_params[1:]);

	for i in range(len(decyear)):
		east_detrended.append(Data0.dE[i]-(east_model[i]) );
		north_detrended.append(Data0.dN[i]-(north_model[i]) );
		vert_detrended.append(Data0.dU[i]-(vert_model[i]) );
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=north_detrended, dE=east_detrended, dU=vert_detrended, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
	return newData;



def pair_gps_model(gps_data, model_data):
	# Takes two time series objects, and returns two paired time series objects. 
	# It could be that GPS has days that model doesn't, or the other way around. 
	dtarray=[]; dE_gps=[]; dN_gps=[]; dU_gps=[]; Se_gps=[]; Sn_gps=[]; Su_gps=[];
	dE_model=[]; dN_model=[]; dU_model=[]; Se_model=[]; Sn_model=[]; Su_model=[]; 
	gps_data=remove_nans(gps_data);
	model_data=remove_nans(model_data);
	for i in range(len(gps_data.dtarray)):
		if gps_data.dtarray[i] in model_data.dtarray:
			idx = model_data.dtarray.index(gps_data.dtarray[i]);  # where is this datetime object in the model array? 
			dtarray.append(gps_data.dtarray[i]);
			dE_gps.append(gps_data.dE[i]);
			dN_gps.append(gps_data.dN[i]);
			dU_gps.append(gps_data.dU[i]);
			Se_gps.append(gps_data.Se[i]);
			Sn_gps.append(gps_data.Sn[i]);
			Su_gps.append(gps_data.Su[i]);
			dE_model.append(model_data.dE[idx]);
			dN_model.append(model_data.dN[idx]);
			dU_model.append(model_data.dU[idx]);
			Se_model.append(model_data.Se[idx]);
			Sn_model.append(model_data.Sn[idx]);
			Su_model.append(model_data.Su[idx]);
	paired_gps = Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=dtarray, dE=dE_gps, dN=dN_gps, dU=dU_gps, Se=Se_gps, Sn=Sn_gps, Su=Su_gps, EQtimes=gps_data.EQtimes);
	paired_model = Timeseries(name=model_data.name, coords=model_data.coords, dtarray=dtarray, dE=dE_model, dN=dN_model, dU=dU_model, Se=Se_model, Sn=Sn_model, Su=Su_model, EQtimes=model_data.EQtimes);
	return [paired_gps, paired_model];


def pair_gps_model_keeping_gps(gps_data, model_data):
	# Takes two time series objects, and returns two time series objects. 
	# It keeps all the data from the first one. 
	# It generates a model_data with length that matches the GPS. 
	dE_model=[]; dN_model=[]; dU_model=[]; Se_model=[]; Sn_model=[]; Su_model=[]; dtarray=[]; 
	gps_data=remove_nans(gps_data);
	model_data=remove_nans(model_data);
	for i in range(len(gps_data.dtarray)):
		if gps_data.dtarray[i] in model_data.dtarray:
			idx = model_data.dtarray.index(gps_data.dtarray[i]);  # where is this datetime object in the model array? 
			dtarray.append(gps_data.dtarray[i]);
			dE_model.append(model_data.dE[idx]);
			dN_model.append(model_data.dN[idx]);
			dU_model.append(model_data.dU[idx]);
			Se_model.append(model_data.Se[idx]);
			Sn_model.append(model_data.Sn[idx]);
			Su_model.append(model_data.Su[idx]);
		else:
			dtarray.append(gps_data.dtarray[i]);  # if we can't find it, then we put filler model. 
			dE_model.append(0);
			dN_model.append(0);
			dU_model.append(0);
			Se_model.append(0);
			Sn_model.append(0);
			Su_model.append(0);
	paired_model = Timeseries(name=model_data.name, coords=model_data.coords, dtarray=dtarray, dE=dE_model, dN=dN_model, dU=dU_model, Se=Se_model, Sn=Sn_model, Su=Su_model, EQtimes=model_data.EQtimes);
	return [gps_data, paired_model];



def get_referenced_data(roving_station_data, base_station_data):
	# Takes a time series object and removes the motion of a base station (another time series object)
	# If there's a starttime, then we will solve for a best-fitting model offset at the starttime. 
	# This is used when subtracting models
	dtarray=[]; dE_gps=[]; dN_gps=[]; dU_gps=[]; Se_gps=[]; Sn_gps=[]; Su_gps=[];
	roving_station_data=remove_nans(roving_station_data);	
	for i in range(len(roving_station_data.dtarray)):
		if roving_station_data.dtarray[i] in base_station_data.dtarray:
			idx = base_station_data.dtarray.index(roving_station_data.dtarray[i]);  # where is this datetime object in the model array? 
			dtarray.append(roving_station_data.dtarray[i]);
			dE_gps.append(roving_station_data.dE[i] - base_station_data.dE[idx]);
			dN_gps.append(roving_station_data.dN[i] - base_station_data.dN[idx]);
			dU_gps.append(roving_station_data.dU[i] - base_station_data.dU[idx]);
			Se_gps.append(roving_station_data.Se[i]);
			Sn_gps.append(roving_station_data.Sn[i]);
			Su_gps.append(roving_station_data.Su[i]);
	gps_relative = Timeseries(name=roving_station_data.name, coords=roving_station_data.coords, dtarray=dtarray, dE=dE_gps, dN=dN_gps, dU=dU_gps, Se=Se_gps, Sn=Sn_gps, Su=Su_gps, EQtimes=roving_station_data.EQtimes);
	return gps_relative; 

def remove_constant(Data0, east_offset, north_offset, vert_offset):
	# Subtract a constant number from each data array in a time series object
	temp_east=[];
	temp_north=[];
	temp_vert=[];
	for i in range(len(Data0.dtarray)):
		temp_east.append(Data0.dE[i]-east_offset);
		temp_north.append(Data0.dN[i]-north_offset);
		temp_vert.append(Data0.dU[i]-vert_offset);
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=temp_north, dE=temp_east, dU=temp_vert, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
	return newData;



# FUTURE FEATURES: 
def rotate_data():
	return;







# -------------------------------------------- # 
# FUNCTIONS THAT TAKE TIME SERIES OBJECTS #
# AND RETURN SCALARS OR VALUES # 
# -------------------------------------------- # 

def get_slope(Data0, starttime=[], endtime=[],missing_fraction=0.6):
	# Model the data with a best-fit y = mx + b. 
	# Returns six numbers: e_slope, n_slope, v_slope, e_std, n_std, v_std

	# Defensive programming
	error_flag, starttime, endtime = basic_defensive_programming(Data0, starttime, endtime);
	if error_flag:
		return [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan];

	# Cut to desired window, and remove nans
	mydtarray=[]; myeast=[]; mynorth=[]; myup=[];

	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime and ~np.isnan(Data0.dE[i]) and ~np.isnan(Data0.dN[i]) and ~np.isnan(Data0.dU[i]):
			mydtarray.append(Data0.dtarray[i]);
			myeast.append(Data0.dE[i]);
			mynorth.append(Data0.dN[i]);
			myup.append(Data0.dU[i]);

	# More defensive programming
	if len(mydtarray)<=2:
		print("ERROR: no time array for station %s. Returning Nan" % Data0.name);
		return [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan];	
	time_duration=mydtarray[-1]-mydtarray[0];
	if time_duration.days<270:
		print("ERROR: using much less than one year of data to estimate parameters for station %s. Returning Nan" % Data0.name);
		return [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan];
	if len(myeast)<time_duration.days*missing_fraction:
		print("ERROR: Most of the data is missing to estimate parameters for station %s. Returning Nan" % Data0.name);
		return [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan];	


	# doing the inversion here, since it's only one line.
	decyear=get_float_times(mydtarray);
	east_coef=np.polyfit(decyear,myeast,1);
	north_coef=np.polyfit(decyear,mynorth,1);
	vert_coef=np.polyfit(decyear,myup,1);
	east_slope=east_coef[0];
	north_slope=north_coef[0];
	vert_slope=vert_coef[0];


	# How bad is the fit to the line? 
	east_trend=[east_coef[0]*x + east_coef[1] for x in decyear];
	east_detrended=[myeast[i]-east_trend[i] for i in range(len(myeast))];
	east_std = np.std(east_detrended);
	north_trend=[north_coef[0]*x + north_coef[1] for x in decyear];
	north_detrended=[mynorth[i]-north_trend[i] for i in range(len(mynorth))];
	north_std = np.std(north_detrended);
	vert_trend=[vert_coef[0]*x + vert_coef[1] for x in decyear];
	vert_detrended=[myup[i]-vert_trend[i] for i in range(len(myup))];
	vert_std = np.std(vert_detrended);	

	return [east_slope, north_slope, vert_slope, east_std, north_std, vert_std];


def get_slope_unc(dataObj, starttime, endtime):
	# Uses the Allan Variance of Rates. 
	dataObj=impose_time_limits(dataObj, starttime, endtime);
	x = get_float_times(dataObj.dtarray);
	y = dataObj.dE;
	sig= dataObj.Se;
	params, covm = lssq_model_errors.AVR(x, y, sig, verbose=0);
	slope = params[0];
	Esigma = np.sqrt(covm[0][0]);
	y = dataObj.dN;
	sig= dataObj.Sn;
	params, covm = lssq_model_errors.AVR(x, y, sig, verbose=0);
	slope = params[0];
	Nsigma = np.sqrt(covm[0][0]);
	y = dataObj.dU;
	sig= dataObj.Su;
	params, covm = lssq_model_errors.AVR(x, y, sig, verbose=0);
	slope = params[0];
	Usigma = np.sqrt(covm[0][0]);
	return [Esigma, Nsigma, Usigma]; 


def get_linear_annual_semiannual(Data0, starttime=[], endtime=[],critical_len=365):
	# The critical_len parameter allows us to manually switch this function for both GPS and GRACE time series in GPS format
	# Model the data with a best-fit GPS = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt) + E*t + F; 

	# Defensive programming
	error_flag, starttime, endtime = basic_defensive_programming(Data0, starttime, endtime);
	if error_flag:
		east_params=[np.nan,0,0,0,0];  north_params=[np.nan,0,0,0,0]; up_params=[np.nan,0,0,0,0];
		return [east_params, north_params, vert_params];
	
	# Cut to desired time window, and remove nans.
	mydtarray=[]; myeast=[]; mynorth=[]; myup=[];
	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime and ~np.isnan(Data0.dE[i]):
			mydtarray.append(Data0.dtarray[i]);
			myeast.append(Data0.dE[i]);
			mynorth.append(Data0.dN[i]);
			myup.append(Data0.dU[i]);

	duration = mydtarray[-1]-mydtarray[0];
	if duration.days<critical_len:
		print("ERROR: using less than one year of data to estimate parameters for station %s. Returning Nan" % Data0.name);
		east_params=[np.nan,0,0,0,0];  north_params=[np.nan,0,0,0,0]; up_params=[np.nan,0,0,0,0];
		return [east_params, north_params, up_params];

	decyear=get_float_times(mydtarray);	
	east_params_unordered=invert_linear_annual_semiannual(decyear,myeast);
	north_params_unordered=invert_linear_annual_semiannual(decyear, mynorth);
	vert_params_unordered=invert_linear_annual_semiannual(decyear, myup);

	# The definition for returning parameters, consistent with Noel's reporting:
	# slope, a2(cos), a1(sin), s2, s1. 
	east_params=[east_params_unordered[4], east_params_unordered[0], east_params_unordered[1], east_params_unordered[2], east_params_unordered[3]];
	north_params=[north_params_unordered[4], north_params_unordered[0], north_params_unordered[1], north_params_unordered[2], north_params_unordered[3]];
	vert_params=[vert_params_unordered[4], vert_params_unordered[0], vert_params_unordered[1], vert_params_unordered[2], vert_params_unordered[3]];

	return [east_params, north_params, vert_params];



def get_means(Data0, starttime=[], endtime=[]):
	# Return the average value of the time series between starttime and endtime
	# Can be used to set offsets for plotting, etc. 

	# Defensive programming
	error_flag, starttime, endtime = basic_defensive_programming(Data0, starttime, endtime);
	if error_flag:
		return [np.nan, np.nan, np.nan];

	# Cut to desired window, and remove nans
	mydtarray=[]; myeast=[]; mynorth=[]; myup=[];
	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime and ~np.isnan(Data0.dE[i]) and ~np.isnan(Data0.dN[i]) and ~np.isnan(Data0.dU[i]):
			mydtarray.append(Data0.dtarray[i]);
			myeast.append(Data0.dE[i]);
			mynorth.append(Data0.dN[i]);
			myup.append(Data0.dU[i]);	

	return [np.nanmean(myeast), np.nanmean(mynorth), np.nanmean(myup)];


def get_logfunction(Data0, eqtime, starttime=[], endtime=[]):
	# y = B + Alog(1+t/tau);
	# Useful for postseismic transients. 
	# Should match the construct-function function. 
	# The t is in decyear
	def func(t, a, b, tau):
		return b + a*np.log(1+t/tau);

	float_times = get_relative_times(Data0.dtarray, dt.datetime.strptime("2010-04-04","%Y-%m-%d"));  # in days
	e_params, ecov = curve_fit(func, float_times, Data0.dE);
	n_params, ecov = curve_fit(func, float_times, Data0.dN);
	u_params, ecov = curve_fit(func, float_times, Data0.dU);
	return [e_params, n_params, u_params];


def get_values_at_date(Data0, selected_date, num_days=10):
	# At a selected date, pull out the east, north, and up values at that date. 
	if selected_date in Data0.dtarray:
		idx = Data0.dtarray.index(selected_date);
		e_value=np.nanmean(Data0.dE[idx:idx+num_days]);
		n_value=np.nanmean(Data0.dN[idx:idx+num_days]);
		u_value=np.nanmean(Data0.dU[idx:idx+num_days]);
	else:
		print("Error: requested date %s not found in dtarray" % dt.datetime.strftime(selected_date, "%Y-%m-%d") );
		[e_value, n_value, u_value] = [np.nan, np.nan, np.nan];
	return e_value, n_value, u_value; 


# -------------------------------------------- # 
# MISCELLANEOUS FUNCTIONS
# SUCH AS MATH AND TIME OPERATIONS
# -------------------------------------------- # 


def basic_defensive_programming(Data0, starttime, endtime):
	# Check for all sorts of nasty things. 
	error_flag = 0;

	if len(Data0.dtarray)==0:
		print("Error: length of dtarray is 0 for station %s. Returning nan" % Data0.name);
		error_flag=1;
		return error_flag, starttime, endtime;

	if starttime==[]:
		starttime=Data0.dtarray[0];
	if endtime==[]:
		endtime=Data0.dtarray[-1];

	starttime_proper=starttime;
	endtime_proper=endtime;

	if starttime<Data0.dtarray[0]:
		starttime_proper=Data0.dtarray[0];
	if endtime>Data0.dtarray[-1]:
		endttime_proper=Data0.dtarray[-1];

	if endtime<Data0.dtarray[0]:
		print("Error: end time before start of array for station %s. Returning Nan" % Data0.name);
		error_flag=1;
	if starttime>Data0.dtarray[-1]:
		print("Error: start time after end of array for station %s. Returning Nan" % Data0.name);
		error_flag=1;

	return error_flag, starttime_proper, endtime_proper; 



def invert_linear_annual_semiannual(decyear,data):
	"""
	Take a time series and fit a best-fitting linear least squares equation: 
	GPS = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt) + E*t + F; 
	Here we also solve for a linear trend as well. 
	"""
	design_matrix=[];
	w = 2*np.pi / 1.0;  
	for t in decyear:
		design_matrix.append([np.cos(w*t), np.sin(w*t), np.cos(2*w*t), np.sin(2*w*t), t, 1]);
	design_matrix= np.array(design_matrix);
	params = np.dot(np.linalg.inv(np.dot(design_matrix.T, design_matrix)), np.dot(design_matrix.T, data));
	return params;


def get_float_times(datetimes):
	floats=[];
	for item in datetimes:
		floats.append(get_float_time(item));
	return floats;

def get_float_time(datetime_item):
	temp=datetime_item.strftime("%Y %j");
	temp=temp.split();
	year = temp[0]
	last_day_of_year = dt.datetime.strptime(year+"1231","%Y%m%d");
	num_days_this_year = float(last_day_of_year.strftime("%Y %j").split()[1]);   #This is either 365 or 366
	floats = (float(temp[0])+float(temp[1])/num_days_this_year);
	return floats;

def get_relative_times(datetimes, origin_dt):
	floats=[];
	for item in datetimes:
		floats.append(get_relative_time(item, origin_dt));
	return floats;

def get_relative_time(datetime_item, origin_dt):
	timedelta = datetime_item - origin_dt;
	relative_time = timedelta.days;
	return relative_time;

def reltime_to_dt(relative_time, origin_dt):
	timedelta = dt.timedelta(days=relative_time);
	my_date = origin_dt+timedelta;
	return my_date;

def float_to_dt(float_time):
	# Example: 2014.194 --> datetime object
	fractional_year=str(1+int(365.24*(float_time-np.floor(float_time))));  # something like 004, 204, 321, etc. 
	if len(fractional_year)==1:
		fractional_year='00'+fractional_year;
	elif len(fractional_year)==2:
		fractional_year='0'+fractional_year;
	if fractional_year=='367' or fractional_year=='366':
		fractional_year='365';
	myyear = str(int(np.floor(float_time)));  # something like 2014
	my_date = dt.datetime.strptime(myyear+fractional_year,"%Y%j");
	return my_date;

def get_daily_dtarray(starttime, endtime):
	# Return a datetime array that spans starttime to endtime with daily intervals
	i=0;
	dtarray=[];
	dtarray.append(starttime);
	while dtarray[-1]<endtime:
		i=i+1;
		dtarray.append(starttime+dt.timedelta(days=i));
	return dtarray;	

def yrnum2datetime(yearnums, starttime):
	# This function will take a set of dates, in decimal years past a certain date, 
	# and convert it into normal datetime objects. 
	# The input and output vectors will be exactly the same length
	dtarray = [];
	for i in range(len(yearnums)):
		myyr = yearnums[i];
		dtarray.append(starttime+dt.timedelta(days=myyr*365.24));
	return dtarray;	

def add_two_unc_quadrature(unc1, unc2):
	return np.sqrt(unc1*unc1 + unc2*unc2);





# -------------------------------------------- # 
# FUNCTIONS THAT TAKE PARAMETERS
# AND RETURN Y=F(X) ARRAYS
# -------------------------------------------- # 


def linear_annual_semiannual_function(decyear, fit_params):
	"""
	Given curve parameters and a set of observation times, build the function y = f(x). 
	Model consists of GPS_V = E*t + Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt);
	"""
	model_def = [];
	w = 2*np.pi / 1.0; 
	for t in decyear:
		model_def.append( fit_params[0]*t + (fit_params[1]*np.cos(w*t)) + (fit_params[2]*np.sin(w*t)) + (fit_params[3]*np.cos(2*w*t)) + (fit_params[4]*np.sin(2*w*t)) );
	return model_def;


def annual_semiannual_only_function(decyear, fit_params):
	"""
	Given curve parameters and a set of observation times, build the function y = f(x). 
	Model consists of GPS_V = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt);
	"""
	model_def = [];
	w = 2*np.pi / 1.0; 
	for t in decyear:
		model_def.append( (fit_params[0]*np.cos(w*t)) + (fit_params[1]*np.sin(w*t)) + (fit_params[2]*np.cos(2*w*t)) + (fit_params[3]*np.sin(2*w*t)) );
	return model_def;


def annual_only_function(decyear, fit_params):
	"""
	Given curve parameters and a set of observation times, build the function y = f(x). 
	Model consists of GPS_V = Acos(wt) + Bsin(wt); 
	"""
	model_def = [];
	w = 2*np.pi / 1.0; 
	for t in decyear:
		model_def.append( (fit_params[0]*np.cos(w*t)) + (fit_params[1]*np.sin(w*t)) );
	return model_def;


def construct_log_function(decday, fit_params):
	# Function has functional form:
	# y = b + a*np.log(1+t/tau);
	# fit params = [a, b, tau];
	# The x axis is the same as the get_log_function()
	a=fit_params[0]; 
	b=fit_params[1];
	tau=fit_params[2];
	model_def = [];
	for i in range(len(decday)):
		model_def.append(b + a*np.log(1+decday[i]/tau));
	return model_def;



