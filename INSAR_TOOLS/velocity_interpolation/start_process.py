"""
The purpose of this script is to compute a least-squares velocity field from the PBO dataset. 

"""

import collections
import glob
import numpy as np 
import datetime as dt 
import matplotlib.pyplot as plt 
import os 
from subprocess import call 

Params=collections.namedtuple("Params",['data_directory','data_center','file_list','output_file'])
data_collection = collections.namedtuple("data_collection",['name','coords','yyyymmdd','dN', 'dE','dU','Sn','Se','Su']);


def make_vel_file():
	myParams = super_configure();
	for filename in myParams.file_list:
		myData = inputs(filename);
		vels   = process(myParams,myData);
		outputs(myParams,myData.name, myData.coords,vels);
	return;


# ----------------- CONFIGURE ------------------- # 
def super_configure():
	data_directory="../pbo_stations/"
	data_center="pbo"  # pbo or unm or cwu
	reference_frame="nam08"
	output_file=reference_frame+"_vel_field.txt";
	if os.path.exists(output_file):
		call(["rm",output_file],shell=False);
	file_list=glob.glob(data_directory+"????."+data_center+".final_"+reference_frame+".pos");
	myParams=Params(data_directory=data_directory, data_center=data_center,file_list=file_list,output_file=output_file);
	return myParams;



# ------------------ INPUTS -------------------- # 
def inputs(filename):
	[yyyymmdd, Nlat, Elong, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, skiprows=37, unpack=True,usecols=(0,12,13,15,16,17,18,19,20));
	specific_file=filename.split('/')[-1];
	myData=data_collection(name=specific_file[0:4],coords=[Elong[0]-360, Nlat[0]], yyyymmdd=yyyymmdd, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su);
	return myData;


# ------------------ PROCESSING -------------------- # 
def process(myParams,myData):
	decyear=get_decyear(myData.yyyymmdd)
	if np.max(decyear) - np.min(decyear) < 3:
		east_velocity=np.NaN;
		north_velocity=np.NaN;
		vertical_velocity=np.NaN;
	else:
		gps_e_params = fit_linear_annual_semiannual(decyear, myData.dE);
		gps_n_params = fit_linear_annual_semiannual(decyear, myData.dN);
		gps_v_params = fit_linear_annual_semiannual(decyear, myData.dU);
		east_velocity= gps_e_params[4]*1000;
		north_velocity= gps_n_params[4]*1000;
		vertical_velocity=gps_v_params[4]*1000;		# Change to millimeters
	return [east_velocity, north_velocity, vertical_velocity];

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
def outputs(myParams, name, coords, vels):
	f1=open(myParams.output_file,'a');
	f1.write("%s %f %f %f %f %f\n" % (name, coords[0], coords[1], vels[0], vels[1], vels[2]) );
	f1.close();
	return;







