# August 2018
# This is a toolbox that operates on Timeseries objects. 
# Its purpose is to deal with antenna offsets and earthquake offsets


import numpy as np 
import collections
import subprocess
import datetime as dt
import sys
from scipy import signal
import gps_io_functions

# A line for referencing the namedtuple definition. 
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm




def parse_antenna_table_pbo(table):
	e_offsets=[]; n_offsets=[]; u_offsets=[]; evdts=[];
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
		e_offsets.append(float(words[8]));  # in m
		n_offsets.append(float(words[6]));
		u_offsets.append(float(words[10]));
		evdts.append(dt.datetime.strptime(yyyy+mm+dd,"%Y%m%d"));	
	return [e_offsets, n_offsets, u_offsets, evdts];


def parse_earthquake_table_pbo(table):
	e_offsets=[]; n_offsets=[]; u_offsets=[]; evdts=[];
	tablesplit=table.split('\n');
	for item in tablesplit:  # for each earthquake
		if len(item)==0:
			continue;  # if we're at the end, move on. 
		words = item.split();
		filename=words[0];
		e_offsets.append(float(words[3]));  # in m
		n_offsets.append(float(words[4]));
		u_offsets.append(float(words[8]));
		evdate = filename.split('/')[-1];
		evdate = evdate[4:10];
		year=evdate[0:2];
		month=evdate[2:4];
		day=evdate[4:6];
		year="20"+year;
		evdts.append(dt.datetime.strptime(year+month+day,"%Y%m%d"));
	return [e_offsets, n_offsets, u_offsets, evdts];


def parse_table_unr(table,offset_type):
	# Here we extract all the antenna or earthquake offsets from the UNR table. 
	# Offset_type is 'antenna'[1] or 'eq'[2]
	evdts=[];
	if offset_type=='antenna':
		desired_code='1';  # antenna offsets
	else:
		desired_code='2';  # earthquakes

	tablesplit=table.split('\n');
	for item in tablesplit:
		if len(item)==0:
			continue;
		words=item.split();
		datestring=words[1];
		evtype=words[2];
		if evtype==desired_code:
			mydt=get_datetime_from_unrfile(datestring);
			evdts.append(mydt);
	return evdts;

def get_datetime_from_unrfile(input_string):
	# Turns something like "12FEB13" into datetime for 2012-02-13
	year=input_string[0:2];
	if int(year)>=80:
		year='19'+year;
	else:
		year='20'+year;
	mydt = dt.datetime.strptime(year+input_string[2:],"%Y%b%d");
	return mydt;



def solve_for_offsets(ts_object, offset_times):
	# Here we solve for all the offsets at a given time, which is necessary for UNR data.
	# Offset_times is a list of datetime objects. 

	index_array=[];
	e_offsets=[]; n_offsets=[]; u_offsets=[]; evdts_new=[];
	number_of_days = 10;

	for i in range(len(offset_times)):  # Find where the object exists in the array
		myindex=np.where(ts_object.dtarray==np.datetime64(offset_times[i]));
		if len(myindex[0])>0:
			# print(ts_object.dtarray[myindex[0][0]]);
			# Find the datetime object that actually matches the time of the offset. 
			index_array.append(myindex[0][0]);  # the i'th measurement is where the offset happens. 
			evdts_new.append(ts_object.dtarray[myindex[0][0]]);

	for i in index_array:
		mean_e_before=np.nanmean(ts_object.dE[i-1-number_of_days:i-1]);
		mean_e_after=np.nanmean(ts_object.dE[i+1:i+1+number_of_days]);
		e_offsets.append(mean_e_after-mean_e_before);

		mean_n_before=np.nanmean(ts_object.dN[i-1-number_of_days:i-1]);
		mean_n_after=np.nanmean(ts_object.dN[i+1:i+1+number_of_days]);
		n_offsets.append(mean_n_after-mean_n_before);

		mean_u_before=np.nanmean(ts_object.dU[i-1-number_of_days:i-1]);
		mean_u_after=np.nanmean(ts_object.dU[i+1:i+1+number_of_days]);
		u_offsets.append(mean_u_after-mean_u_before);

	return [e_offsets, n_offsets, u_offsets, evdts_new];




def remove_offsets_compute(Data0, e_offsets, n_offsets, u_offsets, evdts):
	if len(e_offsets)==0:
		return Data0;
	newdtarray=[]; newdN=[]; newdE=[]; newdU=[];
	# Removing offsets
	for i in range(len(Data0.dtarray)):
		# For each day...
		tempE=Data0.dE[i];
		tempN=Data0.dN[i];
		tempU=Data0.dU[i];
		for j in range(len(evdts)):
			# print("removing %f mm from east at %s" % (e_offset[j], evdt[j]));
			if Data0.dtarray[i]>=evdts[j]:
				tempE=tempE-e_offsets[j];
				tempN=tempN-n_offsets[j];
				tempU=tempU-u_offsets[j];
		newdtarray.append(Data0.dtarray[i]);
		newdE.append(tempE);
		newdN.append(tempN);
		newdU.append(tempU);
	
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=[]);
	return newData;



def remove_antenna_offsets(Data0, offsets_dir, datasource='pbo'):
	# Parse UNR or PBO table, solve for offsets if necessary, remove offsets. 
	station=Data0.name;
	print("Offset table for station %s:" % (station) );	

	# Read the offset table
	if datasource=='pbo':
		try:
			table = subprocess.check_output("grep "+station+" "+offsets_dir+"*.off",shell=True);
		except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
			table=[];
	elif datasource=='unr':
		try:
			table = subprocess.check_output("grep "+station+" "+offsets_dir+"UNR_steps.txt",shell=True);
		except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
			table=[];
	table=table.decode();  # needed when switching to python 3
	print(table);

	# Get the values of offsets. 
	if datasource=='pbo':
		[e_offsets, n_offsets, u_offsets, evdts]=parse_antenna_table_pbo(table);
	elif datasource=='unr':
		evdts=parse_table_unr(table, 'antenna');
		[e_offsets, n_offsets, u_offsets, evdts]=solve_for_offsets(Data0, evdts);

	# Actually remove offsets
	newData=remove_offsets_compute(Data0, e_offsets, n_offsets, u_offsets, evdts);  # the actual subtraction of offsets. 
	return newData; 



def remove_earthquakes(Data0, earthquakes_dir, offsets_dir, datasource='pbo'):
	# Building earthquake table
	station=Data0.name;
	print("Earthquake table for station %s:" % (station) );	

	# Read the offset table
	if datasource=='pbo':
		try:
			table = subprocess.check_output("grep "+station+" "+earthquakes_dir+"*kalts.evt",shell=True);
		except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
			table=[];
	elif datasource=='unr':
		try:
			table = subprocess.check_output("grep "+station+" "+offsets_dir+"UNR_steps.txt",shell=True);
		except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
			table=[];
	table=table.decode(); # needed when switching to python 3
	print(table);

	if datasource=='pbo':
		[e_offsets, n_offsets, u_offsets, evdts]=parse_earthquake_table_pbo(table);
	elif datasource=='unr':
		evdts=parse_table_unr(table,'eq');
		[e_offsets, n_offsets, u_offsets, evdts]=solve_for_offsets(Data0, evdts);  

	newData=remove_offsets_compute(Data0, e_offsets, n_offsets, u_offsets, evdts);  # the actual subtraction of offsets. 
	# Replacing the EQtimes with the proper ones. 
	newData=Timeseries(name=newData.name, coords=newData.coords, dtarray=newData.dtarray, dN=newData.dN, dE=newData.dE, dU=newData.dU, Sn=newData.Sn, Se=newData.Se, Su=newData.Su, EQtimes=evdts);
	return newData;


