# Functions to map, filter, and reduce generic GPS time series


import numpy as np 
import collections
import subprocess
import datetime as dt 
from scipy import signal

# 
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


# -------------------------------------------- # 
# Functions that return modified timeseries
# These all operate on Timeseries objects. 
# -------------------------------------------- # 

def remove_earthquakes(Data0, earthquakes_dir):
	# Building earthquake table
	station=Data0.name;
	table = subprocess.check_output("grep "+station+" "+earthquakes_dir+"*kalts.evt",shell=True);
	print("building earthquake table for station %s ..." % (station) );
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

	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
	return newData;


def detrend_data(Data0):
	east_detrended=signal.detrend(Data0.dE);
	north_detrended=signal.detrend(Data0.dN);
	vert_detrended=signal.detrend(Data0.dU);
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=north_detrended, dE=east_detrended, dU=vert_detrended, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);	
	return newData;

def rotate_data():
	return;



# -------------------------------------------- # 
# FUNCTIONS THAT RETURN SCALARS OR VALUES # 
# -------------------------------------------- # 



