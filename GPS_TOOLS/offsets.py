# August 2018
# This is a toolbox that operates on Timeseries objects. 
# Its purpose is to deal with antenna offsets and earthquake offsets


import numpy as np 
import collections
import datetime as dt

# A line for referencing the namedtuple definition. 
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Offsets = collections.namedtuple("Offsets",['e_offsets', 'n_offsets', 'u_offsets', 'evdts']);

def remove_offsets(Data0, offsets_obj):
	# the actual subtraction of offsets. 
	if offsets_obj==[]:
		return Data0;
	if len(offsets_obj.e_offsets)==0:
		return Data0;
	newdtarray=[]; newdN=[]; newdE=[]; newdU=[];

	# Removing offsets
	for i in range(len(Data0.dtarray)):
		# For each day...
		tempE=Data0.dE[i];
		tempN=Data0.dN[i];
		tempU=Data0.dU[i];
		for j in range(len(offsets_obj.evdts)):
			# print("removing %f mm from east at %s" % (offsets_obj.e_offsets[j], offsets_obj.evdts[j]));
			if Data0.dtarray[i]==offsets_obj.evdts[j]:   # removing the date of the offset directly (it can be messed up)
				tempE=np.nan;
				tempN=np.nan;
				tempU=np.nan;
			if Data0.dtarray[i]>offsets_obj.evdts[j]:
				tempE=tempE-offsets_obj.e_offsets[j];
				tempN=tempN-offsets_obj.n_offsets[j];
				tempU=tempU-offsets_obj.u_offsets[j];
		newdtarray.append(Data0.dtarray[i]);
		newdE.append(tempE);
		newdN.append(tempN);
		newdU.append(tempU);
	
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
	return newData;



def fit_offset(dtarray, data, interval, offset_num_days):
	# Loop through the array and find dates that are used for offset calculation. 
	# This is done for one component, like east
	# The offset time can be a day or an interval (day is just repeated twice.)
	before_indeces = [];
	after_indeces = [];

	# Find the indeces of nearby days
	for i in range(len(dtarray)):
		deltat_start=dtarray[i]-interval[0];  # the beginning of the interval
		deltat_end = dtarray[i]-interval[1];  # the end of the interval
		if deltat_start.days >= -offset_num_days and deltat_start.days<=0: 
			before_indeces.append(i);
		if deltat_end.days <= offset_num_days and deltat_end.days >= 0: 
			after_indeces.append(i);

	# Identify the value of the offset. 
	if before_indeces==[] or after_indeces==[] or len(before_indeces)==1 or len(after_indeces)==1:
		offset=0;
		print("Warning: no data before or after offset at %s. Returning offset=0" % (dt.datetime.strftime(interval[0],"%Y-%m-%d")) );
	else:
		before_mean= np.nanmean( [data[x] for x in before_indeces] );
		after_mean = np.nanmean( [data[x] for x in after_indeces] );
		offset=after_mean-before_mean;
		if offset==np.nan:
			print("Warning: np.nan offset found. Returning 0");
			offset=0;
	return offset;


def solve_for_offsets(ts_object, offset_times):
	# Here we solve for all the offsets at a given time, which is necessary for UNR data.
	# Offset_times is a list of datetime objects with unique dates. 
	index_array=[];
	e_offsets=[]; n_offsets=[]; u_offsets=[]; evdts_new=[];
	number_of_days = 10;

	for i in range(len(offset_times)):
		e_offset = fit_offset(ts_object.dtarray, ts_object.dE, [offset_times[i], offset_times[i]], number_of_days);
		n_offset = fit_offset(ts_object.dtarray, ts_object.dN, [offset_times[i], offset_times[i]], number_of_days);
		u_offset = fit_offset(ts_object.dtarray, ts_object.dU, [offset_times[i], offset_times[i]], number_of_days);
		e_offsets.append(e_offset);
		n_offsets.append(n_offset);
		u_offsets.append(u_offset);
		evdts_new.append(offset_times[i]);

	Offset_obj = Offsets(e_offsets=e_offsets, n_offsets=n_offsets, u_offsets=u_offsets, evdts=evdts_new);
	return Offset_obj;


def get_empty_offsets():
	Offset = Offsets(e_offsets=[], n_offsets=[], u_offsets=[], evdts=[]);
	return Offset;

def print_offset_object(Offset_obj):
	print("Total offset object:");
	for i in range(len(Offset_obj.e_offsets)):
		print("%s: %.4f mmE, %.4f mmN, %.4f mmU" % (dt.datetime.strftime(Offset_obj.evdts[i],"%Y-%m-%d"), Offset_obj.e_offsets[i], Offset_obj.n_offsets[i], Offset_obj.u_offsets[i]) );
	return;


