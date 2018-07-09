"""
The purpose of this script is to view the vertical velocity field from the PBO dataset. 

"""

import collections
import glob
import numpy as np 
import datetime as dt 
import matplotlib.pyplot as plt 
import os 
from subprocess import call 


Params=collections.namedtuple("Params",['data_directory','data_center','file_list','lonlimits','latlimits','output_file']);
data_collection = collections.namedtuple("data_collection",['name','coords','yyyymmdd','dN', 'dE','dU','Sn','Se','Su']);


def compute_linear_fits():
	myParams = configure();
	myData = inputs(myParams);
	[evel,nvel,vvel] = compute(myParams,myData);
	outputs(myParams,myData,evel,nvel,vvel);
	return;


# ----------------- CONFIGURE ------------------- # 
def configure():
	data_directory="../pbo_stations/"
	data_center="pbo"  # pbo or unm or cwu
	output_file="vel_field.txt";
	file_list=glob.glob(data_directory+"????."+data_center+".final_igs08.pos");
	lonlimits=[-125, -119];
	latlimits=[38, 43.3];
	myParams=Params(data_directory=data_directory, data_center=data_center,file_list=file_list,lonlimits=lonlimits, latlimits=latlimits,output_file=output_file);
	return myParams;



# ------------------ INPUTS -------------------- # 
def inputs(myParams):  # inputs for all stations. 
	names=[]; yrstrings=[]; coords=[]; dNs=[]; dEs=[]; dUs=[]; Sns=[]; Ses=[]; Sus=[];

	for filename in myParams.file_list:
		[yyyymmdd, Nlat, Elong, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, skiprows=37, unpack=True,usecols=(0,12,13,15,16,17,18,19,20));
		name=filename.split('/')[-1];
		names.append(name);
		yrstrings.append(yyyymmdd);
		coords.append([Elong[0]-360, Nlat[0]]);
		dNs.append(dN);
		dEs.append(dE);
		dUs.append(dU);
		Sns.append(Sn);
		Ses.append(Se);
		Sus.append(Su);
	
	myData=data_collection(name=names,coords=coords, yyyymmdd=yrstrings, dN=dNs, dE=dEs, dU=dUs, Sn=Sns, Se=Ses, Su=Sus);
	return myData;


# ------------------ PROCESSING -------------------- # 
def compute(myParams,myData):
	evel=[]; nvel=[]; vvel=[];

	for i in range(len(myData.name)):
		decyear=get_decyear(myData.yyyymmdd[i])
		if np.max(decyear) - np.min(decyear) < 3:
			station_vvel=np.NaN;
			station_evel=np.NaN;
			station_nvel=np.NaN;
		else:
			gps_e_params = fit_linear_annual_semiannual(decyear, myData.dE[i]);
			station_evel=gps_e_params[4]*1000;		# Change to millimeters
			gps_n_params = fit_linear_annual_semiannual(decyear, myData.dN[i]);
			station_nvel=gps_n_params[4]*1000;		# Change to millimeters
			gps_v_params = fit_linear_annual_semiannual(decyear, myData.dU[i]);
			station_vvel=gps_v_params[4]*1000;		# Change to millimeters	
		evel.append(station_evel);
		nvel.append(station_nvel);
		vvel.append(station_vvel);
	return [evel, nvel, vvel];

def get_decyear(yyyymmdd):
	decyear=[];
	for t in yyyymmdd:
		mydt = dt.datetime.strptime(str(int(t)), '%Y%m%d');
		year=mydt.strftime('%Y');
		jday=mydt.strftime('%j');
		decyear.append(float(year)+float(jday)/365.24);
	return decyear;

def fit_linear_annual_semiannual(decyear,data):
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



# ------------------ OUTPUTS -------------------- # 
def outputs(myParams,myData,evel,nvel,vvel):
	f1=open(myParams.output_file,'a');
	for i in range(len(evel)):
		f1.write("%s %f %f %f %f %f \n" % (myData.name[i], myData.coords[i][0], myData.coords[i][1], evel[i], nvel[i], vvel[i]) );
	f1.close();
	plotting_outputs(myParams);
	return;


def plotting_outputs(myParams):
	[lon, lat, vvel] = np.loadtxt(myParams.output_file,usecols=(1,2,3),unpack=True);
	[ca_lon, ca_lat] = np.loadtxt('california_bdr',unpack=True);
	[nv_lon, nv_lat] = np.loadtxt('nevada_bdr',unpack=True);
	[or_lon, or_lat] = np.loadtxt('oregon_bdr',unpack=True);

	# Set upper and lower bounds on uplift values
	Vbounds = [-3, 3];	
	for i in range(len(vvel)):
		if vvel[i]<Vbounds[0]:
			vvel[i]=Vbounds[0];
		if vvel[i]>Vbounds[1]:
			vvel[i]=Vbounds[1];	

	plt.figure();
	plt.plot(ca_lon,ca_lat,'k');
	plt.plot(nv_lon,nv_lat,'k');
	plt.plot(or_lon,or_lat,'k');
	plt.scatter(lon,lat,s=40, c=vvel,edgecolor=None)
	cb = plt.colorbar();
	cb.set_label('Velocity (mm/yr)');

	plt.xlim(myParams.lonlimits);
	plt.ylim(myParams.latlimits);
	plt.savefig("vvel_field_ccal.jpg");
	
	plt.close();
	return;




