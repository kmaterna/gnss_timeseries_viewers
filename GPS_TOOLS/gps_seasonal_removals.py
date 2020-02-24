# A set of functions that take TimeSeries objects and return other TimeSeries objects
# Without seasonal components. 

import numpy as np 
import matplotlib.pyplot as plt
import collections, sys
import subprocess
import datetime as dt 
import glob
import gps_ts_functions
import notch_filter
import grace_ts_functions
import gps_io_functions


Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


# Make a detrended/modeled version of this time series. 
def make_detrended_ts(Data, seasonals_remove, seasonals_type, remove_trend=1, 
	fit_table="../../GPS_POS_DATA/Velocity_Files/Bartlow_interETSvels.txt",
	grace_dir="../../GPS_POS_DATA/GRACE_loading_model/",
	STL_dir="../../GPS_POS_DATA/STL_models/",
	gldas_dir="../../GPS_POS_DATA/PBO_Hydro/GLDAS/",
	nldas_dir="../../GPS_POS_DATA/PBO_Hydro/NLDAS/", 
	lakes_dir="../../GPS_POS_DATA/Lake_loading/",
	lsdm_dir="../../GPS_POS_DATA/PBO_Hydro/LSDM/"):
	# Once we have removed earthquake steps... 
	# The purpose of this function is to generate a version of the time series that has been detrended and optionally seasonal-removed, 
	# Where the seasonal fitting (if necessary) and detrending happen in the same function. 
	# There are options for the seasonal removal (least squares, notch filter, grace, etc.)
	# Fit params definition: slope, a2(cos), a1(sin), s2, s1. 

	# Here we are asking to invert the data for linear and seasonal components
	if seasonals_remove==0:
		print("Not removing seasonals.");
		east_params=[0,0,0,0,0];  north_params=[0,0,0,0,0]; up_params=[0,0,0,0,0];
		[east_vel, north_vel, up_vel,esig,nsig,usig]=gps_ts_functions.get_slope(Data);
		east_params[0]=east_vel; north_params[0]=north_vel; up_params[0]=up_vel;
		trend_out=gps_ts_functions.detrend_data_by_value(Data, east_params, north_params, up_params);
		trend_in = Data;

	else:  # Going into different forms of seasonal removal. 
		print("Removing seasonals by %s method." % seasonals_type);
		if seasonals_type=='lssq':
			[east_params, north_params, up_params]=gps_ts_functions.get_linear_annual_semiannual(Data);
			trend_out= gps_ts_functions.detrend_data_by_value(Data, east_params, north_params, up_params);
			trend_in = gps_ts_functions.remove_seasonal_by_value(Data, east_params, north_params, up_params);

		elif seasonals_type=='noel':
			[east_params, north_params, up_params] = look_up_seasonal_coefs(Data.name, fit_table);
			trend_out= gps_ts_functions.detrend_data_by_value(Data, east_params, north_params, up_params);	
			trend_in = Data; # Broken because I don't really use this option anymore. 
		
		elif seasonals_type=='notch':
			trend_out, trend_in = remove_seasonals_by_notch(Data);

		elif seasonals_type=='grace':
			trend_out, trend_in = remove_seasonals_by_GRACE(Data,grace_dir);

		elif seasonals_type=='stl':
			trend_out, trend_in = remove_seasonals_by_STL(Data, STL_dir);

		elif seasonals_type=='nldas':
			trend_out, trend_in = remove_seasonals_by_hydro(Data, nldas_dir);

		elif seasonals_type=='nldas_scaled':
			trend_out, trend_in = remove_seasonals_by_hydro(Data, nldas_dir, scaling=True);

		elif seasonals_type=='gldas':
			trend_out, trend_in = remove_seasonals_by_hydro(Data, gldas_dir);

		elif seasonals_type=='lsdm':
			trend_out, trend_in = remove_seasonals_by_german_load(Data, lsdm_dir);

		elif seasonals_type=='oroville':
			trend_out, trend_in = remove_seasonals_by_lakes(Data, lakes_dir, 'oroville');

		elif seasonals_type=='shasta':
			trend_out, trend_in = remove_seasonals_by_lakes(Data, lakes_dir, 'shasta');

		else:
			print("Error: %s not supported as a seasonal removal type" % seasonals_type);
			print("The supported types are: lssq, noel, grace, notch, and stl");
			print("Exiting!\n");
			sys.exit(1);

	if remove_trend==0:
		return trend_in;
	else:
		return trend_out;



# This function is the function that actually gets called by other programs. 
# It operates on a TimeSeries object and returns another TimeSeries object
def remove_seasonals_by_notch(Data):
	# Using Sang-Ho's notch filter script to remove power at frequencies corresponding to 1 year and 6 months. 
	# We are also removing a linear trend in this step. 

	Data=gps_ts_functions.remove_nans(Data);

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

	decyear=gps_ts_functions.get_float_times(Data.dtarray);
	dE_detrended=np.zeros(np.shape(Data.dE)); dN_detrended=np.zeros(np.shape(Data.dN)); dU_detrended=np.zeros(np.shape(Data.dU));
	dE_trended=np.zeros(np.shape(Data.dE)); dN_trended=np.zeros(np.shape(Data.dN)); dU_trended=np.zeros(np.shape(Data.dU));

	# East
	x = Data.dE;
	dE_filt = notch_filter.notchfilt(x,fs,fn1,Bn1,filtfiltopt=True);
	dE_filt = notch_filter.notchfilt(dE_filt,fs,fn2,Bn2,filtfiltopt=True);
	east_coef=np.polyfit(decyear,dE_filt,1)[0];
	for i in range(len(dE_filt)):
		dE_detrended[i]=dE_filt[i]-east_coef*decyear[i] - (dE_filt[0]-east_coef*decyear[0]);
		dE_trended[i]  =dE_filt[i]; 

	# North
	x = Data.dN;
	dN_filt = notch_filter.notchfilt(x,fs,fn1,Bn1,filtfiltopt=True);
	dN_filt = notch_filter.notchfilt(dN_filt,fs,fn2,Bn2,filtfiltopt=True);
	north_coef=np.polyfit(decyear,dN_filt,1)[0];
	for i in range(len(dN_filt)):
		dN_detrended[i]=dN_filt[i]-north_coef*decyear[i] - (dN_filt[0]-north_coef*decyear[0]);
		dN_trended[i]  =dN_filt[i];

	# Up
	x = Data.dU;
	dU_filt = notch_filter.notchfilt(x,fs,fn1,Bn1,filtfiltopt=True);
	dU_filt = notch_filter.notchfilt(dU_filt,fs,fn2,Bn2,filtfiltopt=True);
	vert_coef=np.polyfit(decyear,dU_filt,1)[0];
	for i in range(len(dU_filt)):
		dU_detrended[i]=dU_filt[i]-vert_coef*decyear[i] - (dU_filt[0]-vert_coef*decyear[0]);
		dU_trended[i] = dU_filt[i];

	detrended=Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=dN_detrended, dE=dE_detrended, dU=dU_detrended, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
	trended  =Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=dN_trended, dE=dE_trended, dU=dU_trended, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
	return detrended, trended;


def remove_seasonals_by_STL(Data, STL_dir):
	# Has an issue: Not sure if it returns trended data. 
	# Right now only returns detrended data. 

	# First check if a pre-computed file exists

	filename=STL_dir+Data.name+"_STL_30.txt";
	recompute=0;
	try:
		ifile=open(filename);
	except FileNotFoundError:
		print("Warning! STL not found for %s" % Data.name);
		recompute=1;

	if recompute==0:
		# If a precomputed file exists... 
		[dE, dN, dU, Se, Sn, Su]=np.loadtxt(filename,unpack=True,usecols=(1,2,3,4,5,6));
		final_dtarray=[];
		for line in ifile:
			dtstring=line.split()[0];
			final_dtarray.append(dt.datetime.strptime(dtstring,"%Y%m%d"));
		Data=Timeseries(name=Data.name, coords=Data.coords, dtarray=final_dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su, EQtimes=Data.EQtimes);

	# ELSE: WE NEED TO RECOMPUTE
	else:
		print("We did not find a pre-computed array, so we are re-computing STL. ");

		# Preprocess data: remove nans, fill in gaps. 
		Data=gps_ts_functions.remove_nans(Data);
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
		decyear=gps_ts_functions.get_float_times(new_dtarray);		

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
		
		# Write the file so that we don't recompute it next time. 
		output_stl(Data,STL_dir);

	return Data, Data;

def output_stl(Data, outdir):
	ofile=open(outdir+Data.name+"_STL_30.txt",'w');
	for i in range(len(Data.dtarray)):
		timestamp=dt.datetime.strftime(Data.dtarray[i],"%Y%m%d");
		ofile.write("%s %f %f %f %f %f %f\n" % (timestamp, Data.dE[i], Data.dN[i], Data.dU[i], Data.Se[i], Data.Sn[i], Data.Su[i]) );
	ofile.close();
	return;


def preprocess_stl(dtarray, data_column, uncertainties):
	# fill in gaps, and make to a multiple of 365 day cycles. 

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


def remove_seasonals_by_hydro(Data, hydro_dir, scaling=False):
	station=Data.name;
	files=glob.glob(hydro_dir+station.lower()+'*.hyd');
	# filename=hydro_dir+station.lower()+'_noah10_gldas2.hyd';
	if len(files)>0:
		filename=files[0];
	else:
		filename='';

	try:
		ifile=open(filename);
		print("Opening file %s" % filename);	
	except FileNotFoundError:
		print("Error! Hydro file not found for %s" % Data.name);
		placeholder = np.full_like(Data.dtarray, np.nan, dtype=np.double)
		wimpyObj=Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=placeholder, dE=placeholder, dU=placeholder, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
		print("returning placeholder object");
		return wimpyObj, wimpyObj;  # 1 = error code. 

	# Read the hydro model and pair it to the GPS
	[hydro_data]=gps_io_functions.read_pbo_hydro_file(filename);
	
	# Clean up and pair data
	Data=gps_ts_functions.remove_nans(Data);
	# hydro_data=gps_ts_functions.remove_nans(hydro_data);  # this may or may not be necessary. 
	[gps_data, hydro_data] = gps_ts_functions.pair_gps_model(Data, hydro_data);  # matched in terms of dtarray. 

	if scaling==True:
		[east_gps, north_gps, vert_gps] = gps_ts_functions.get_linear_annual_semiannual(gps_data);
		[east_hydro, north_hydro, vert_hydro] = gps_ts_functions.get_linear_annual_semiannual(hydro_data);
		gps_amp = np.sqrt(vert_gps[1]*vert_gps[1] + vert_gps[2]*vert_gps[2]);
		hydro_amp = np.sqrt(vert_hydro[1]*vert_hydro[1] + vert_hydro[2]*vert_hydro[2]);
		if hydro_amp == 0.0:
			print("ERROR! NLDAS amplitude is exactly 0!!  You should probably fix this. ");
			placeholder = np.full_like(Data.dtarray, np.nan, dtype=np.double)
			wimpyObj=Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=placeholder, dE=placeholder, dU=placeholder, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
			print("returning placeholder object");
			return wimpyObj, wimpyObj;  # 1 = error code. 			
		scale_factor=gps_amp/hydro_amp;
		# print("GPS Amplitude is %.2f mm" % gps_amp);
		# print("NLDAS Amplitude is %.2f mm" % hydro_amp);
		print("NLDAS scaling factor is %.2f" % scale_factor);
	else:
		scale_factor=1;

	#  Subtract the model from the data. 
	dE_filt=[]; dN_filt=[]; dU_filt=[];
	for i in range(len(gps_data.dtarray)):
		dE_filt.append(gps_data.dE[i]-scale_factor*hydro_data.dE[i]);
		dN_filt.append(gps_data.dN[i]-scale_factor*hydro_data.dN[i]);
		dU_filt.append(gps_data.dU[i]-scale_factor*hydro_data.dU[i]);

	# A Simple detrending
	decyear = gps_ts_functions.get_float_times(gps_data.dtarray);
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

	corrected_object=Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_detrended, dN=dN_detrended, dU=dU_detrended, Se=gps_data.Se, Sn=gps_data.Sn, Su=gps_data.Su, EQtimes=gps_data.EQtimes);
	trended=Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_filt, dN=dN_filt, dU=dU_filt, Se=gps_data.Se, Sn=gps_data.Sn, Su=gps_data.Su, EQtimes=gps_data.EQtimes);
	return corrected_object, trended;



def remove_seasonals_by_german_load(Data, lsdm_dir):
	station=Data.name;
	files=glob.glob(lsdm_dir+station+'*.txt');

	if len(files)>0:
		filename=files[0];
	else:
		filename='';

	try:
		ifile=open(filename);
		print("Opening file %s" % filename);	
	except FileNotFoundError:
		print("Error! LSDM file not found for %s" % Data.name);
		placeholder = np.full_like(Data.dtarray, np.nan, dtype=np.double)
		wimpyObj=Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=placeholder, dE=placeholder, dU=placeholder, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
		print("returning placeholder object");
		return wimpyObj, wimpyObj;  # 1 = error code. 

	# Read the hydro model and pair it to the GPS
	[hydro_data]=gps_io_functions.read_lsdm_file(filename);
	
	# Clean up and pair data
	Data=gps_ts_functions.remove_nans(Data);
	hydro_data=gps_ts_functions.remove_nans(hydro_data);  # this may or may not be necessary. 
	[gps_data, hydro_data] = gps_ts_functions.pair_gps_model(Data, hydro_data);  # matched in terms of dtarray. 

	#  Subtract the model from the data. 
	dE_filt=[]; dN_filt=[]; dU_filt=[];
	for i in range(len(gps_data.dtarray)):
		dE_filt.append(gps_data.dE[i]-hydro_data.dE[i]);
		dN_filt.append(gps_data.dN[i]-hydro_data.dN[i]);
		dU_filt.append(gps_data.dU[i]-hydro_data.dU[i]);

	# A Simple detrending
	decyear = gps_ts_functions.get_float_times(gps_data.dtarray);
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

	detrended=Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_detrended, dN=dN_detrended, dU=dU_detrended, Se=gps_data.Se, Sn=gps_data.Sn, Su=gps_data.Su, EQtimes=gps_data.EQtimes);
	trended=Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_filt, dN=dN_filt, dU=dU_filt, Se=gps_data.Se, Sn=gps_data.Sn, Su=gps_data.Su, EQtimes=gps_data.EQtimes);
	return detrended, trended;


def remove_seasonals_by_lakes(Data, lakes_dir, lake_name):
	station=Data.name;
	files=glob.glob(lakes_dir+station+"_"+lake_name+"*.txt");
	if len(files)==1:
		print("Finding lake deformation data in %s" % files[0]);
	else:
		print("Error! Lake %s file not found for %s" % (lake_name, Data.name));
		placeholder = np.full_like(Data.dtarray, np.nan, dtype=np.double)
		wimpyObj=Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=placeholder, dE=placeholder, dU=placeholder, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes);
		print("returning placeholder object");
		return wimpyObj, wimpyObj;  # 1 = error code. 

	loading_ts = read_loading_ts(files[0]);
	[GPS_paired, loading_paired] = gps_ts_functions.pair_gps_model(Data, loading_ts);
	dE=[GPS_paired.dE[i]-loading_paired.dE[i] for i in range(len(GPS_paired.dE))];
	dN=[GPS_paired.dN[i]-loading_paired.dN[i] for i in range(len(GPS_paired.dN))];
	dU=[GPS_paired.dU[i]-loading_paired.dU[i] for i in range(len(GPS_paired.dU))];

	decyear=gps_ts_functions.get_float_times(GPS_paired.dtarray);
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

	corrected_object = Timeseries(name=Data.name, coords=Data.coords, dtarray=GPS_paired.dtarray, dE=dE_detrended, dN=dN_detrended, dU=dU_detrended, Se=GPS_paired.Se, Sn=GPS_paired.Sn, Su=GPS_paired.Su, EQtimes=Data.EQtimes);
	trended = Timeseries(name=Data.name, coords=Data.coords, dtarray=GPS_paired.dtarray, dE=dE, dN=dN, dU=dU, Se=GPS_paired.Se, Sn=GPS_paired.Sn, Su=GPS_paired.Su, EQtimes=Data.EQtimes);
	return corrected_object, trended;


def read_loading_ts(infile):
	dtarray = []; u = []; v = []; w = []; 
	ifile=open(infile,'r');
	for line in ifile:
		temp=line.split();
		dtarray.append(dt.datetime.strptime(temp[0],"%Y-%m-%d"));
		u.append(float(temp[4]));
		v.append(float(temp[5]));
		w.append(float(temp[6]));
	S=np.zeros(np.shape(u));
	# HERE WE WILL MAKE A NEW DATA OBJECT
	ifile.close();
	loading_defo = Timeseries(name='', coords=[], dtarray=dtarray, dE=u, dN=v, dU=w, Sn=S, Se=S, Su=S, EQtimes=[]);
	return loading_defo;


def look_up_seasonal_coefs(name,table_file):
	[E, N, U, Ea1, Na1, Ua1, Ea2, Na2, Ua2, Es1, Ns1, Us1, Es2, Ns2, Us2]=gps_io_functions.read_noel_file_station(table_file,name);
	east_params=[E, Ea2, Ea1, Es2, Es1];
	north_params=[N, Na2, Na1, Ns2, Na2];
	up_params=[U, Ua2, Ua1, Us2, Us1];
	return [east_params, north_params, up_params];


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
		print("returning placeholder object");
		return wimpyObj;  # 1 = error code. 

	# If the station has been pre-computed with GRACE:
	Data=gps_ts_functions.remove_nans(Data);
	grace_model=grace_ts_functions.input_GRACE_individual_station(grace_dir+"scaled_"+Data.name+"_PREM_model_ts.txt");
	my_paired_ts = grace_ts_functions.pair_GPSGRACE(Data, grace_model);
	decyear = gps_ts_functions.get_float_times(my_paired_ts.dtarray);

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

	detrended=Timeseries(name=Data.name, coords=Data.coords, dtarray=my_paired_ts.dtarray, dN=dN_detrended, dE=dE_detrended, dU=dU_detrended, Sn=my_paired_ts.N_err, Se=my_paired_ts.E_err, Su=my_paired_ts.V_err, EQtimes=Data.EQtimes);
	trended=Timeseries(name=Data.name, coords=Data.coords, dtarray=my_paired_ts.dtarray, dN=dN_filt, dE=dE_filt, dU=dU_filt, Sn=my_paired_ts.N_err, Se=my_paired_ts.E_err, Su=my_paired_ts.V_err, EQtimes=Data.EQtimes);
	return detrended, trended; # 0 = successful completion



