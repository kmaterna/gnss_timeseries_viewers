# Read the inputs for PBO and UNR files

import numpy as np 
import collections
import subprocess, sys, os
import datetime as dt
import gps_io_functions
import offsets

Offsets = collections.namedtuple("Offsets",['e_offsets', 'n_offsets', 'u_offsets', 'evdts']);
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


# ----------------------------------------
# DRIVERS, CONFIGURE, AND FILE MASHING ---
# ----------------------------------------

def get_unr(station):
	unr_filename="../../GPS_POS_DATA/UNR_Data/"+station+".NA12.tenv3"
	unr_coords="../../GPS_POS_DATA/UNR_DATA/UNR_coords_july2018.txt"
	offsets_dir="../../GPS_POS_DATA/Offsets/"
	[myData]=gps_io_functions.read_UNR_magnet_file(unr_filename, unr_coords);  # UNR data format
	Offsets = get_unr_offsets(myData, station, offsets_dir);
	Earthquakes = get_unr_earthquakes(myData, station, offsets_dir);
	return [myData, Offsets, Earthquakes];


def get_pbo(station):
	pbo_filename="../../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos"
	pbo_earthquakes_dir="../../GPS_POS_DATA/PBO_Event_Files/"
	offsets_dir="../../GPS_POS_DATA/Offsets/"
	[myData]=gps_io_functions.read_pbo_pos_file(pbo_filename);  # PBO data format
	Offsets = get_pbo_offsets(station, offsets_dir);
	Earthquakes = get_pbo_earthquakes(station, pbo_earthquakes_dir);
	return [myData, Offsets, Earthquakes];


# Based on whether a file exists in certain directories or not, 
# Return the 'pbo' or 'unr' datasource that we should be using. 
def determine_datasource(station, input_datasource):
	unr_filename="../../GPS_POS_DATA/UNR_Data/"+station+".NA12.tenv3";
	pbo_filename="../../GPS_POS_DATA/PBO_Data/"+station+".pbo.final_nam08.pos";
	if input_datasource=='pbo' and os.path.isfile(pbo_filename):
		print("Using PBO file as input data. ");
		datasource='pbo'
	elif input_datasource=='pbo' and os.path.isfile(unr_filename):
		print("Using UNR as input data because PBO data file doesn't exist");
		datasource='unr';
	elif input_datasource=='unr' and not os.path.isfile(unr_filename):
		print("Error! Cannot find file in UNR database.");
		sys.exit(1);
	elif input_datasource != 'unr' and input_datasource != 'pbo':
		print("Error! Invalid input datasource");
		sys.exit(1);
	else:
		print("Cannot find input file for station %s ; exiting..." % station);
		sys.exit(1);
	return datasource;


# ----------------------------------------------- 
# THE GUTS --------------------------------------
# -----------------------------------------------

def get_unr_offsets(Data0, station, offsets_dir):
	print("Offset table for station %s:" % (station) );	
	try:
		table = subprocess.check_output("grep "+station+" "+offsets_dir+"UNR_steps.txt",shell=True);
	except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
		table=[];
	table=table.decode();  # needed when switching to python 3
	print(table);
	evdts=parse_table_unr(table, 'antenna');
	[e_offsets, n_offsets, u_offsets, evdts]=solve_for_offsets(Data0, evdts);
	UNR_offsets=Offsets(e_offsets=e_offsets, n_offsets=n_offsets, u_offsets=u_offsets, evdts=evdts);
	return UNR_offsets;


def get_unr_earthquakes(Data0, station, offsets_dir):
	print("Earthquake table for station %s:" % (station) );	
	try:
		table = subprocess.check_output("grep "+station+" "+offsets_dir+"UNR_steps.txt",shell=True);
	except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
		table=[];
	table=table.decode(); # needed when switching to python 3
	print(table);
	evdts=parse_table_unr(table,'eq');
	[e_offsets, n_offsets, u_offsets, evdts]=solve_for_offsets(Data0, evdts); 
	UNR_earthquakes=Offsets(e_offsets=e_offsets, n_offsets=n_offsets, u_offsets=u_offsets, evdts=evdts);
	return UNR_earthquakes;



def get_pbo_offsets(station, offsets_dir):
	print("Offset table for station %s:" % (station) );	
	try:
		table = subprocess.check_output("grep "+station+" "+offsets_dir+"*.off",shell=True);
	except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
		table=[];
	table=table.decode();  # needed when switching to python 3
	print(table);
	[e_offsets, n_offsets, u_offsets, evdts]=parse_antenna_table_pbo(table);
	PBO_offsets=Offsets(e_offsets=e_offsets, n_offsets=n_offsets, u_offsets=u_offsets, evdts=evdts);
	return PBO_offsets;

def get_pbo_earthquakes(station, earthquakes_dir):
	print("Earthquake table for station %s:" % (station) );	
	# Read the offset table
	try:
		table = subprocess.check_output("grep "+station+" "+earthquakes_dir+"*kalts.evt",shell=True);
	except subprocess.CalledProcessError as grepexc:  # if we have no earthquakes in the event files... 
		table=[];
	table=table.decode(); # needed when switching to python 3
	print(table);
	[e_offsets, n_offsets, u_offsets, evdts]=parse_earthquake_table_pbo(table);
	PBO_earthquakes=Offsets(e_offsets=e_offsets, n_offsets=n_offsets, u_offsets=u_offsets, evdts=evdts);
	return PBO_earthquakes;


#
# TABLE INPUTS --------------------------- 
# 

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

#
# TABLE COMPUTE --------------------------- 
# 

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

