# May/June 2018
# This is a toolbox that operates on Timeseries objects. 
# Contains functions to map, filter, and reduce generic GPS time series


import numpy as np 
import collections
import subprocess
import datetime as dt 
import sys
from scipy import signal

# A line for referencing the namedtuple definition. 
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


# -------------------------------------------- # 
# FUNCTIONS THAT TAKE TIME SERIES OBJECTS # 
# AND RETURN OTHER TIME SERIES OBJECTS # 
# -------------------------------------------- # 

def remove_offsets(Data0, offsets_dir):
	station=Data0.name;
	e_offset=[]; n_offset=[]; u_offset=[]; evdt=[];
	print("Offset table for station %s:" % (station) );	

	try:
		table = subprocess.check_output("grep "+station+" "+offsets_dir+"*.off",shell=True);
	except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
		table=[];

	if len(table)==0:
		return Data0;

	table_rows=table.split('\n');
	for line in table_rows:
		if "EQ" in line:
			continue;
		else:
			print(line);
		if len(line)==0:
			continue;  # if we're at the end, move on. 
		words = line.split();
		site=words[0];
		yyyy=words[1];
		mm=words[2];
		dd=words[3];
		e_offset.append(float(words[8]));  # in m
		n_offset.append(float(words[6]));
		u_offset.append(float(words[10]));
		evdt.append(dt.datetime.strptime(yyyy+mm+dd,"%Y%m%d"));

	newdtarray=[]; newdN=[]; newdE=[]; newdU=[];
	# Removing offsets
	for i in range(len(Data0.dtarray)):
		# For each day...
		tempE=Data0.dE[i];
		tempN=Data0.dN[i];
		tempU=Data0.dU[i];
		for j in range(len(evdt)):
			# print("removing %f mm from east at %s" % (e_offset[j], evdt[j]));
			if Data0.dtarray[i]>=evdt[j]:
				tempE=tempE-e_offset[j];
				tempN=tempN-n_offset[j];
				tempU=tempU-u_offset[j];
		newdtarray.append(Data0.dtarray[i]);
		newdE.append(tempE);
		newdN.append(tempN);
		newdU.append(tempU);
	
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=evdt);

	return newData; 


def remove_earthquakes(Data0, earthquakes_dir):
	# Building earthquake table
	station=Data0.name;

	try:
		table = subprocess.check_output("grep "+station+" "+earthquakes_dir+"*kalts.evt",shell=True);
	except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
		table=[];

	if len(table)==0:
		return Data0;

	print("Earthquake table for station %s:" % (station) );
	print(table);
	e_offset=[]; n_offset=[]; u_offset=[]; evdt=[];
	tablesplit=table.split('\n');
	for item in tablesplit:  # for each earthquake
		if len(item)==0:
			continue;  # if we're at the end, move on. 
		words = item.split();
		filename=words[0];
		e_offset.append(float(words[3]));  # in m
		n_offset.append(float(words[4]));
		u_offset.append(float(words[8]));
		evdate = filename.split('/')[-1];
		evdate = evdate[4:10];
		year=evdate[0:2];
		month=evdate[2:4];
		day=evdate[4:6];
		year="20"+year;
		evdt.append(dt.datetime.strptime(year+month+day,"%Y%m%d"));

	newdtarray=[]; newdN=[]; newdE=[]; newdU=[];
	# Removing offsets
	for i in range(len(Data0.dtarray)):
		# For each day...
		tempE=Data0.dE[i];
		tempN=Data0.dN[i];
		tempU=Data0.dU[i];
		for j in range(len(evdt)):
			# print("removing %f mm from east at %s" % (e_offset[j], evdt[j]));
			if Data0.dtarray[i]>=evdt[j]:
				tempE=tempE-e_offset[j];
				tempN=tempN-n_offset[j];
				tempU=tempU-u_offset[j];
		newdtarray.append(Data0.dtarray[i]);
		newdE.append(tempE);
		newdN.append(tempN);
		newdU.append(tempU);
	
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=evdt);
	return newData;


def remove_outliers(Data0, outliers_def):
	medfilt_e=signal.medfilt(Data0.dE, 35);
	medfilt_n=signal.medfilt(Data0.dN, 35);
	medfilt_u=signal.medfilt(Data0.dU, 35);

	newdtarray=[]; newdN=[]; newdE=[]; newdU=[];
	for i in range(len(medfilt_e)):
		if abs(Data0.dE[i]-medfilt_e[i])<outliers_def and abs(Data0.dN[i]-medfilt_n[i])<outliers_def and abs(Data0.dU[i]-medfilt_u[i])<outliers_def*2:
			newdtarray.append(Data0.dtarray[i]);
			newdE.append(Data0.dE[i]);
			newdN.append(Data0.dN[i]);
			newdU.append(Data0.dU[i]);
		else:
			newdtarray.append(Data0.dtarray[i]);
			newdE.append(np.nan);
			newdN.append(np.nan);
			newdU.append(np.nan);
	
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
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


def detrend_data(Data0):

	# Operates with nans. 
	idx=np.isnan(Data0.dE);
	if(sum(idx))>0:  # if there are nans, please pull them out. 
		Data0=remove_nans(Data0);
	
	decyear=get_float_times(Data0.dtarray);
	east_coef=np.polyfit(decyear,Data0.dE,1);
	north_coef=np.polyfit(decyear,Data0.dN,1);
	vert_coef=np.polyfit(decyear,Data0.dU,1);

	east_detrended=[]; north_detrended=[]; vert_detrended=[];
	for i in range(len(decyear)):
		east_detrended.append(Data0.dE[i]-(east_coef[0]*decyear[i] + east_coef[1]));
		north_detrended.append(Data0.dN[i]-(north_coef[0]*decyear[i] + north_coef[1]));
		vert_detrended.append(Data0.dU[i]-(vert_coef[0]*decyear[i] + vert_coef[1]));
	# print(east_coef[0]);  # these are the slopes
	# print(north_coef[0]);
	# print(vert_coef[0]);
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=north_detrended, dE=east_detrended, dU=vert_detrended, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
	return newData;



def remove_nans(Data0):
	idx=np.isnan(Data0.dE);
	temp_dates=[];
	temp_east=[];
	temp_north=[];
	temp_vert=[];
	temp_Sn=[];
	temp_Se=[];
	temp_Su=[];
	for i in range(len(Data0.dtarray)):
		if idx[i]!=1:
			temp_dates.append(Data0.dtarray[i]);
			temp_east.append(Data0.dE[i]);
			temp_north.append(Data0.dN[i]);
			temp_vert.append(Data0.dU[i]);
			temp_Se.append(Data0.Se[i]);
			temp_Sn.append(Data0.Sn[i]);
			temp_Su.append(Data0.Su[i]);
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=temp_dates, dN=temp_north, dE=temp_east, dU=temp_vert, Sn=temp_Sn, Se=temp_Se, Su=temp_Su, EQtimes=Data0.EQtimes);
	return newData;


def remove_annual_semiannual(Data0):
	decyear=get_float_times(Data0.dtarray);
	[east_params,north_params,up_params]=get_linear_annual_semiannual(Data0);
	east_function = annual_semiannual_only_function(decyear, east_params);
	north_function = annual_semiannual_only_function(decyear, north_params);
	up_function = annual_semiannual_only_function(decyear, up_params);

	temp_east=[]; temp_north=[]; temp_up=[];
	for i in range(len(decyear)):
		temp_east.append(Data0.dE[i]-east_function[i]);
		temp_north.append(Data0.dN[i]-north_function[i]);
		temp_up.append(Data0.dU[i]-up_function[i]);

	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=temp_north, dE=temp_east, dU=temp_up, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes); 
	return newData;



def rotate_data():
	return;





# -------------------------------------------- # 
# FUNCTIONS THAT TAKE TIME SERIES OBJECTS #
# AND RETURN SCALARS OR VALUES # 
# -------------------------------------------- # 

def get_slope(Data0, starttime=[], endtime=[]):
	# Model the data with a best-fit y = mx + b. 
	if starttime==[]:
		starttime=Data0.dtarray[0];
	if endtime==[]:
		endtime=Data0.dtarray[-1];

	# Defensive programming
	if starttime<Data0.dtarray[0]:
		starttime=Data0.dtarray[0];
	if endtime>Data0.dtarray[-1]:
		endttime=Data0.dtarray[-1];
	if endtime<Data0.dtarray[0]:
		print("Error: end time before start of array. Returning Nan");
		return [np.nan,np.nan,np.nan];
	if starttime>Data0.dtarray[-1]:
		print("Error: start time after end of array. Returning Nan");
		return [np.nan,np.nan,np.nan];


	mydtarray=[]; myeast=[]; mynorth=[]; myup=[];
	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime and ~np.isnan(Data0.dE[i]):
			mydtarray.append(Data0.dtarray[i]);
			myeast.append(Data0.dE[i]);
			mynorth.append(Data0.dN[i]);
			myup.append(Data0.dU[i]);

	if len(mydtarray)<365:
		print("Error: using less than one year of data to estimate parameters. Returning 0");
		return [np.nan,np.nan,np.nan];

	decyear=get_float_times(mydtarray);
	east_coef=np.polyfit(decyear,myeast,1);
	north_coef=np.polyfit(decyear,mynorth,1);
	vert_coef=np.polyfit(decyear,myup,1);
	east_slope=east_coef[0];
	north_slope=north_coef[0];
	vert_slope=vert_coef[0];
	return [east_slope, north_slope, vert_slope];


def get_linear_annual_semiannual(Data0, starttime=[], endtime=[]):
	# Model the data with a best-fit GPS = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt) + E*t + F; 
	if starttime==[]:
		starttime=Data0.dtarray[0];
	if endtime==[]:
		endtime=Data0.dtarray[-1];

	# Defensive programming
	if starttime<Data0.dtarray[0]:
		starttime=Data0.dtarray[0];
	if endtime>Data0.dtarray[-1]:
		endttime=Data0.dtarray[-1];
	if endtime<Data0.dtarray[0]:
		print("Error: end time before start of array. Returning 0");
		return [0,0,0];
	if starttime>Data0.dtarray[-1]:
		print("Error: start time after end of array. Returning 0");
		return [0,0,0];

	mydtarray=[]; myeast=[]; mynorth=[]; myup=[];
	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime and ~np.isnan(Data0.dE[i]):
			mydtarray.append(Data0.dtarray[i]);
			myeast.append(Data0.dE[i]);
			mynorth.append(Data0.dN[i]);
			myup.append(Data0.dU[i]);

	if len(mydtarray)<365:
		print("Error: using less than one year of data to estimate parameters. Returning 0");
		return [0,0,0];

	decyear=get_float_times(mydtarray);
	east_params=invert_linear_annual_semiannual(decyear,myeast);
	north_params=invert_linear_annual_semiannual(decyear, mynorth);
	vert_params=invert_linear_annual_semiannual(decyear, myup);
	return [east_params, north_params, vert_params];



def get_float_times(datetimes):
	floats=[];
	for item in datetimes:
		temp=item.strftime("%Y %j");
		temp=temp.split();
		floats.append(float(temp[0])+float(temp[1])/365.24);
	return floats;


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




# -------------------------------------------- # 
# FUNCTIONS THAT TAKE PARAMETERS
# AND RETURN Y=F(X) ARRAYS
# -------------------------------------------- # 


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






