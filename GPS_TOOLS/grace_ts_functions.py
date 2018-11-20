# GRACE FUNCTIONS
import numpy as np
import matplotlib.pyplot as plt 
import collections
import datetime as dt 

Timeseries = collections.namedtuple("Timeseries",[
	'name','coords','dtarray',
	'dN', 'dE','dU',
	'Sn','Se','Su','EQtimes']);  # in mm
GRACE_TS_Array=collections.namedtuple('GRACE_TS',[
	'decyear','dtarray',
	'u','v','w'])
Paired_TS=collections.namedtuple('Paired_TS',[
	'dtarray',
	'north','east','vert',
	'N_err','E_err','V_err',
	'u','v','w']);


def input_GRACE_individual_station(filename):
	# THE GRACE DATA
	try:
		[u, v, w] = np.loadtxt(filename,usecols=range(4,7),unpack=True);
	except FileNotFoundError:
		print("ERROR! Cannot find GRACE model for file %s" % filename);

	u=np.array(u); v=np.array(v); w=np.array(w);
	grace_t = get_grace_datetimes(filename);
	grace_t = [i+dt.timedelta(days=15) for i in grace_t];  # we add 15 days to plot the GRACE data at the center of the bin. 
	grace_decyear=get_float_times(grace_t);  # the decimal years of all the grace obs points. 
	myGRACE_TS=GRACE_TS_Array(decyear=grace_decyear, dtarray=grace_t, u=u, v=v, w=w);
	return myGRACE_TS;


def pair_GPSGRACE(GPS_TS, GRACE_TS):
	# This resamples the GRACE data to match GPS that is within the range of GRACE, and forms a common time axis. 
	gps_decyear=get_float_times(GPS_TS.dtarray)
	decyear=[]; dt=[]; north_gps=[]; east_gps=[]; vert_gps=[]; N_err=[]; E_err=[]; V_err=[]; u=[]; v=[]; w=[];
	for i in range(len(GPS_TS.dtarray)): # this if-statement is happening because GPS is more current than GRACE
		if GPS_TS.dtarray[i]>min(GRACE_TS.dtarray) and GPS_TS.dtarray[i]<max(GRACE_TS.dtarray):
			decyear.append(gps_decyear[i]);
			dt.append(GPS_TS.dtarray[i])
			north_gps.append(GPS_TS.dN[i]);
			east_gps.append(GPS_TS.dE[i]);
			vert_gps.append(GPS_TS.dU[i]);
			N_err.append(GPS_TS.Sn[i]); 
			E_err.append(GPS_TS.Se[i]); 
			V_err.append(GPS_TS.Su[i]); 
	grace_u=np.interp(decyear,GRACE_TS.decyear,GRACE_TS.u);
	grace_v=np.interp(decyear,GRACE_TS.decyear,GRACE_TS.v);
	grace_w=np.interp(decyear,GRACE_TS.decyear,GRACE_TS.w);
	my_paired_ts=Paired_TS(dtarray=dt,north=north_gps,east=east_gps,vert=vert_gps,N_err=N_err,E_err=E_err,V_err=V_err,u=grace_u,v=grace_v,w=grace_w);
	return my_paired_ts;


def get_grace_datetimes(tsfile):
	ifile=open(tsfile,'r');
	dateobjects=[];
	for line in ifile:
		temp=line.split();
		raw_string=temp[0];   # This is in the format 01-Jan-2012_31-Jan-2012
		datestring=raw_string[0:11];
		myobject=dt.datetime.strptime(datestring,'%d-%b-%Y');
		dateobjects.append(myobject);
	ifile.close();
	return dateobjects;

def get_float_times(datetimes):
	floats=[];
	for item in datetimes:
		temp=item.strftime("%Y %j");
		temp=temp.split();
		floats.append(float(temp[0])+float(temp[1])/366.0);
	return floats;


def plot_grace(station_name, filename, out_dir):

	grace_ts = input_GRACE_individual_station(filename);

	plt.figure();
	plt.plot_date(grace_ts.dtarray,grace_ts.u,'-b');
	plt.plot_date(grace_ts.dtarray,grace_ts.v,'-g');
	plt.plot_date(grace_ts.dtarray,grace_ts.w,'-r');
	plt.legend(['east','north','vertical']);
	plt.grid(True);
	plt.xlabel('Time');
	plt.ylabel('Displacement (mm)');
	plt.savefig(out_dir+station_name+"_gracets.eps");

	return;


