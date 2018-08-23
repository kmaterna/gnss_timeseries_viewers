# August 2018
# This is a toolbox that operates on Timeseries objects. 
# Its purpose is to deal with antenna offsets and earthquake offsets


import numpy as np 
import collections
import datetime as dt

# A line for referencing the namedtuple definition. 
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm
Offsets = collections.namedtuple("Offsets",['e_offsets', 'n_offsets', 'u_offsets', 'evdts']);


def remove_antenna_offsets(Data0, offsets_obj):
	# Actually remove offsets
	newData=remove_offsets_compute(Data0, offsets_obj.e_offsets, offsets_obj.n_offsets, offsets_obj.u_offsets, offsets_obj.evdts);  # the actual subtraction of offsets. 
	return newData; 

def remove_earthquakes(Data0, offsets_obj):
	newData=remove_offsets_compute(Data0, offsets_obj.e_offsets, offsets_obj.n_offsets, offsets_obj.u_offsets, offsets_obj.evdts);  # the actual subtraction of offsets. 
	# Replacing the EQtimes with the proper ones. 
	newData=Timeseries(name=newData.name, coords=newData.coords, dtarray=newData.dtarray, dN=newData.dN, dE=newData.dE, dU=newData.dU, Sn=newData.Sn, Se=newData.Se, Su=newData.Su, EQtimes=offsets_obj.evdts);
	return newData;

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



