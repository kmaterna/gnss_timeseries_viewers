# May/June 2018
# This is a toolbox that operates on Timeseries objects. 
# Contains functions to map, filter, and reduce generic GPS time series


import numpy as np 
import collections
import subprocess
import datetime as dt 
import sys
from scipy import signal
import gps_io_functions
import notch_filter
import grace_ts_functions

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


# Make a detrended/modeled version of this time series. 
def make_detrended_option(Data, seasonals_remove, seasonals_type, fit_table="../../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt", grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"):
	# Once we have removed earthquake steps... 
	# The purpose of this function is to generate a version of the time series that has been detrended and optionally seasonal-removed, 
	# Where the seasonal fitting (if necessary) and detrending happen in the same function. 
	# There are options for the seasonal removal (least squares, notch filter, grace, etc.)
	# Fit params definition: slope, a2(cos), a1(sin), s2, s1. 

	east_params=[0,0,0,0,0];  north_params=[0,0,0,0,0]; up_params=[0,0,0,0,0];

	# Here we are asking to invert the data for linear and seasonal components
	if seasonals_type=='fit':
		if seasonals_remove==1:
			[east_params, north_params, up_params]=get_linear_annual_semiannual(Data);
		else:
			[east_vel, north_vel, up_vel]=get_slope(Data);
			east_params[0]=east_vel; north_params[0]=north_vel; up_params[0]=up_vel;
		trend_out=detrend_data_by_value(Data, east_params, north_params, up_params);

	# Here we are removing Noel's fits to the data
	elif seasonals_type=='noel':
		[east_params, north_params, up_params] = look_up_seasonal_coefs(Data, fit_table);
		if seasonals_remove==0:
			east_params[1:-1]=0; north_params[1:-1]=0; up_params[1:-1]=0;  # do not model the seasonal terms
		trend_out=detrend_data_by_value(Data, east_params, north_params, up_params);	
	
	# Here we use notch filters 
	elif seasonals_type=='notch':
		trend_out=remove_seasonals_by_notch(Data);

	# Here we use a pre-computed GRACE loading time series
	elif seasonals_type=='grace':
		trend_out=remove_seasonals_by_GRACE(Data,grace_dir);

	# Here we are doing something else. 
	elif seasonals_type=='stl':
		trend_out=remove_seasonals_by_STL(Data);

	else:
		print("Error: %s not supported as a seasonal removal type" % seasonals_type);


	return trend_out;



def look_up_seasonal_coefs(Data0,table_file):
	[E, N, U, Ea1, Na1, Ua1, Ea2, Na2, Ua2, Es1, Ns1, Us1, Es2, Ns2, Us2]=gps_io_functions.read_noel_file_station(table_file,Data0.name);
	east_params=[E, Ea2, Ea1, Es2, Es1];
	north_params=[N, Na2, Na1, Ns2, Na2];
	up_params=[U, Ua2, Ua1, Us2, Us1];
	return [east_params, north_params, up_params];


def detrend_data_by_value(Data0,east_params,north_params,vert_params):
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


def remove_seasonals_by_notch(Data):
	# Using Sang-Ho's notch filter script to remove power at frequencies corresponding to 1 year and 6 months. 
	# We are also removing a linear trend in this step. 

	Data=remove_nans(Data);

	# Parameters
	# %   x       1-D signal array
	# %   fs      sampling frequency, Hz
	# %   fn      notch frequency, Hz
	# %   Bn      notch bandwidth, Hz
	dt = 1.0;  # one day
	fs = 1/dt;
	fn1 = 1.0/365.24;  # fn = notch frequency, annual
	Bn1 = 0.1*fn1;
	fn2 = 2.0/365.24;  # fn = notch frequency, semiannual
	Bn2 = 0.1*fn2;  # a choice: 10% seems to work well.

	decyear=get_float_times(Data.dtarray);
	dE_detrended=np.zeros(np.shape(Data.dE)); dN_detrended=np.zeros(np.shape(Data.dN)); dU_detrended=np.zeros(np.shape(Data.dU));

	# East
	x = Data.dE;
	dE_filt = notch_filter.notchfilt(x,fs,fn1,Bn1,filtfiltopt=True);
	dE_filt = notch_filter.notchfilt(dE_filt,fs,fn2,Bn2,filtfiltopt=True);
	east_coef=np.polyfit(decyear,dE_filt,1)[0];
	for i in range(len(dE_filt)):
		dE_detrended[i]=dE_filt[i]-east_coef*decyear[i] - (dE_filt[0]-east_coef*decyear[0]);

	# North
	x = Data.dN;
	dN_filt = notch_filter.notchfilt(x,fs,fn1,Bn1,filtfiltopt=True);
	dN_filt = notch_filter.notchfilt(dN_filt,fs,fn2,Bn2,filtfiltopt=True);
	north_coef=np.polyfit(decyear,dN_filt,1)[0];
	for i in range(len(dN_filt)):
		dN_detrended[i]=(dN_filt[i]-north_coef*decyear[i]) - (dN_filt[0]-north_coef*decyear[0]);

	# Up
	x = Data.dU;
	dU_filt = notch_filter.notchfilt(x,fs,fn1,Bn1,filtfiltopt=True);
	dU_filt = notch_filter.notchfilt(dU_filt,fs,fn2,Bn2,filtfiltopt=True);
	vert_coef=np.polyfit(decyear,dU_filt,1)[0];
	for i in range(len(dU_filt)):
		dU_detrended[i]=(dU_filt[i]-vert_coef*decyear[i]) - (dU_filt[0]-vert_coef*decyear[0]);

	newData=Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=dN_detrended, dE=dE_detrended, dU=dU_detrended, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
	return newData;


def remove_seasonals_by_STL(Data):

	# Preprocess data
	[new_dtarray, dE, Se] = preprocess_stl(Data.dtarray, Data.dE, Data.Se);
	[new_dtarray, dN, Sn] = preprocess_stl(Data.dtarray, Data.dN, Data.Sn);
	[new_dtarray, dU, Su] = preprocess_stl(Data.dtarray, Data.dU, Data.Su);

	# Write E, N, U
	ofile=open('raw_ts_data.txt','w');
	for i in range(len(dE)):
		mystring=dt.datetime.strftime(new_dtarray[i],"%Y%m%d");
		ofile.write('%s %f %f %f\n' % (mystring, dE[i], dN[i], dU[i]) );
	ofile.close();

	# Call driver in matlab (read, STL, write)
	subprocess.call(['matlab','-nodisplay','-nosplash','-r','stl_driver'],shell=False);

	# Read / Detrended data
	[dE, dN, dU]=np.loadtxt('filtered_ts_data.txt',unpack=True,usecols=(1,2,3));

	# East, North, Up Detrending
	decyear=get_float_times(new_dtarray);
	dE_detrended=np.zeros(np.shape(dE)); dN_detrended=np.zeros(np.shape(dN)); dU_detrended=np.zeros(np.shape(dU));
	east_coef=np.polyfit(decyear,dE,1)[0];
	for i in range(len(dE)):
		dE_detrended[i]=dE[i]-east_coef*decyear[i] - (dE[0]-east_coef*decyear[0]);
	north_coef=np.polyfit(decyear,dN,1)[0];
	for i in range(len(dN)):
		dN_detrended[i]=(dN[i]-north_coef*decyear[i]) - (dN[0]-north_coef*decyear[0]);
	vert_coef=np.polyfit(decyear,dU,1)[0];
	for i in range(len(dU)):
		dU_detrended[i]=(dU[i]-vert_coef*decyear[i]) - (dU[0]-vert_coef*decyear[0]);

	# Put the gaps back in:
	final_dtarray=[]; final_dE=[]; final_dN=[]; final_dU=[]; final_Se=[]; final_Sn=[]; final_Su=[];
	for i in range(len(new_dtarray)):
		if new_dtarray[i] in Data.dtarray:
			final_dtarray.append(new_dtarray[i]);
			final_dE.append(dE_detrended[i]);
			final_dN.append(dN_detrended[i]);
			final_dU.append(dU_detrended[i]);
			final_Se.append(Se[i]);
			final_Sn.append(Sn[i]);
			final_Su.append(Su[i]);
	final_dE=np.array(final_dE);
	final_dN=np.array(final_dN);
	final_dU=np.array(final_dU);
	final_Se=np.array(final_Se);
	final_Sn=np.array(final_Sn);
	final_Su=np.array(final_Su);

	# Return data
	Data=Timeseries(name=Data.name, coords=Data.coords, dtarray=final_dtarray, dN=final_dN, dE=final_dE, dU=final_dU, Sn=final_Sn, Se=final_Se, Su=final_Su, EQtimes=Data.EQtimes);
	return Data;


def preprocess_stl(dtarray, data_column, uncertainties):
	# Remove nan's, fill in gaps, and make to a multiple of 365 day cycles. 

	new_data_column=[]; new_sig=[]; new_dtarray=[];
	new_data_column.append(data_column[0]);
	new_sig.append(uncertainties[0]);
	new_dtarray.append(dtarray[0]);
	start_date = dtarray[0];
	end_date=dtarray[-1];

	for i in range(1,len(dtarray)):

		start_counter=new_dtarray[-1];  # this is the date where we start. 
		destination_counter=dtarray[i];  # the next datapoint we are going to append until in this loop.

		while start_counter+dt.timedelta(days=1) <= destination_counter:

			# If your next element is consecutive with the old element.
			if start_counter+dt.timedelta(days=1) == destination_counter:	
				new_dtarray.append(dtarray[i]);
				if data_column[i]==np.nan:
					new_data_column.append(new_data_column[-1]);
					new_sig.append(uncertainties[-1]);
				else:
					new_data_column.append(data_column[i]);
					new_sig.append(uncertainties[i]);
				break;

			# If your next element is not consecutive with the old element
			else:
				new_dtarray.append(start_counter+dt.timedelta(days=1));
				new_data_column.append(new_data_column[-1]);
				new_sig.append(uncertainties[-1]);
				start_counter=start_counter+dt.timedelta(days=1);

	# Make the length an integer number of years
	while np.mod(len(new_dtarray),365)>0.01:
		new_dtarray.append(new_dtarray[-1]+dt.timedelta(days=1));
		new_data_column.append(new_data_column[-1]);
		new_sig.append(uncertainties[-1]);

	return [new_dtarray, new_data_column, new_sig];



# Note: 
# Paired_TS=collections.namedtuple('Paired_TS',[
# 	'dtarray',
# 	'north','east','vert',
# 	'N_err','E_err','V_err',
# 	'u','v','w']);

def remove_seasonals_by_GRACE(Data, grace_dir):
	# Here we use pre-computed GRACE load model time series to correct the GPS time series. 
	# We recognize that the horizontals will be bad, and that the resolution of GRACE is coarse.  
	# For these reasons, this is not an important part of the analysis. 
	# Read and interpolate GRACE loading model
	# Subtract the GRACE model
	# Remove a trend from the GPS data
	# Return the object. 

	filename=grace_dir+"scaled_"+Data.name+"_PREM_model_ts.txt";
	try:
		ifile=open(filename);
	except FileNotFoundError:
		print("Error! GRACE not found for %s" % Data.name);
		placeholder = np.full_like(Data.dtarray, np.nan, dtype=np.double)
		wimpyObj=Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=placeholder, dE=placeholder, dU=placeholder, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
		return wimpyObj;

	# If the station has been pre-computed with GRACE:
	Data=remove_nans(Data);
	grace_model=grace_ts_functions.input_GRACE_individual_station(grace_dir+"scaled_"+Data.name+"_PREM_model_ts.txt");
	my_paired_ts = grace_ts_functions.pair_GPSGRACE(Data, grace_model);
	decyear = get_float_times(my_paired_ts.dtarray);

	# Subtract the GRACE object
	dE_filt=[]; dN_filt=[]; dU_filt=[];
	for i in range(len(my_paired_ts.dtarray)):
		dE_filt.append(my_paired_ts.east[i]-my_paired_ts.u[i]);
		dN_filt.append(my_paired_ts.north[i]-my_paired_ts.v[i]);
		dU_filt.append(my_paired_ts.vert[i]-my_paired_ts.w[i]);

	# A Simple detrending
	dE_detrended=np.zeros(np.shape(decyear)); dN_detrended=np.zeros(np.shape(decyear)); dU_detrended=np.zeros(np.shape(decyear));
	east_coef=np.polyfit(decyear,dE_filt,1)[0];
	for i in range(len(dE_filt)):
		dE_detrended[i]=dE_filt[i]-east_coef*decyear[i] - (dE_filt[0]-east_coef*decyear[0]);	
	north_coef=np.polyfit(decyear,dN_filt,1)[0];
	for i in range(len(dN_filt)):
		dN_detrended[i]=(dN_filt[i]-north_coef*decyear[i]) - (dN_filt[0]-north_coef*decyear[0]);
	vert_coef=np.polyfit(decyear,dU_filt,1)[0];
	for i in range(len(dU_filt)):
		dU_detrended[i]=(dU_filt[i]-vert_coef*decyear[i]) - (dU_filt[0]-vert_coef*decyear[0]);

	newData=Timeseries(name=Data.name, coords=Data.coords, dtarray=my_paired_ts.dtarray, dN=dN_detrended, dE=dE_detrended, dU=dU_detrended, Sn=my_paired_ts.N_err, Se=my_paired_ts.E_err, Su=my_paired_ts.V_err, EQtimes=Data.EQtimes);
	return newData;







# FUTURE FEATURES: 
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
		print("Error: end time before start of array for station %s. Returning Nan" % Data0.name);
		return [np.nan,np.nan,np.nan];
	if starttime>Data0.dtarray[-1]:
		print("Error: start time after end of array for station %s. Returning Nan" % Data0.name);
		return [np.nan,np.nan,np.nan];

	# Cut to desired window, and remove nans
	mydtarray=[]; myeast=[]; mynorth=[]; myup=[];
	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime and ~np.isnan(Data0.dE[i]):
			mydtarray.append(Data0.dtarray[i]);
			myeast.append(Data0.dE[i]);
			mynorth.append(Data0.dN[i]);
			myup.append(Data0.dU[i]);

	time_duration=mydtarray[-1]-mydtarray[0];
	if time_duration.days<365:
		print("Error: using less than one year of data to estimate parameters for station %s. Returning 0" % Data0.name);
		return [np.nan,np.nan,np.nan];

	# doing the inversion here, since it's only one line.
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
		print("Error: end time before start of array for station %s. Returning 0" % Data0.name);
		return [0,0,0];
	if starttime>Data0.dtarray[-1]:
		print("Error: start time after end of array for station %s. Returning 0" % Data0.name);
		return [0,0,0];

	# Cut to desired time window, and remove nans.
	mydtarray=[]; myeast=[]; mynorth=[]; myup=[];
	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>=starttime and Data0.dtarray[i]<=endtime and ~np.isnan(Data0.dE[i]):
			mydtarray.append(Data0.dtarray[i]);
			myeast.append(Data0.dE[i]);
			mynorth.append(Data0.dN[i]);
			myup.append(Data0.dU[i]);

	if len(mydtarray)<365:
		print("Error: using less than one year of data to estimate parameters for station %s. Returning 0" % Data0.name);
		return [0,0,0];

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
		temp=item.strftime("%Y %j");
		temp=temp.split();
		floats.append(float(temp[0])+float(temp[1])/365.24);
	return floats;

def get_float_time(datetime_item):
	temp=datetime_item.strftime("%Y %j");
	temp=temp.split();
	floats = (float(temp[0])+float(temp[1])/365.24);
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






