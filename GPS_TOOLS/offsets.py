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
			if Data0.dtarray[i]>=offsets_obj.evdts[j]:
				tempE=tempE-offsets_obj.e_offsets[j];
				tempN=tempN-offsets_obj.n_offsets[j];
				tempU=tempU-offsets_obj.u_offsets[j];
		newdtarray.append(Data0.dtarray[i]);
		newdE.append(tempE);
		newdN.append(tempN);
		newdU.append(tempU);
	
	newData=Timeseries(name=Data0.name, coords=Data0.coords, dtarray=newdtarray, dN=newdN, dE=newdE, dU=newdU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
	return newData;


