"""
The purpose of this script is to view the vertical velocity field from the PBO dataset. 

"""

import collections
import glob
import numpy as np 
import datetime as dt 
import matplotlib.pyplot as plt 
import os 

Params=collections.namedtuple("Params",['data_directory','data_center','file_list','lonlimits','latlimits','output_dir'])
data_collection = collections.namedtuple("data_collection",['name','coords','yyyymmdd','dtarray','decyeararray','dN', 'dE','dU','Sn','Se','Su']);


def stacked_plot():
	myParams = configure();
	myData = inputs(myParams);
	outputs(myParams,myData);
	return;


# ----------------- CONFIGURE ------------------- # 
def configure():
	data_directory="../pbo_stations/"
	data_center="pbo"  # pbo or unm or cwu
	output_dir="";
	file_list=glob.glob(data_directory+"????."+data_center+".final_igs08.pos");
	lonlimits=[-125, -123.7];
	latlimits=[47, 49];
	myParams=Params(data_directory=data_directory, data_center=data_center,file_list=file_list,lonlimits=lonlimits, latlimits=latlimits,output_dir=output_dir);
	return myParams;



# ------------------ INPUTS -------------------- # 
def inputs(myParams):  # inputs for all stations. 
	names=[]; yrstrings=[]; dtarrays=[]; decyeararrays=[]; coords=[]; dNs=[]; dEs=[]; dUs=[]; Sns=[]; Ses=[]; Sus=[];

	for filename in myParams.file_list:
		[yyyymmdd, Nlat, Elong, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, skiprows=37, unpack=True,usecols=(0,12,13,15,16,17,18,19,20));
		if Nlat[0]>myParams.latlimits[0] and Nlat[0]<myParams.latlimits[1]:
			if Elong[0]-360>myParams.lonlimits[0] and Elong[0]-360<myParams.lonlimits[1]:
				
				name=filename.split('/')[-1];
				names.append(name);
				yrstrings.append(yyyymmdd);
				dtarray=[]; decyeararray=[];
				for item in yyyymmdd:
					dtarray.append(dt.datetime.strptime(str(int(item)), "%Y%m%d"))
					decyeararray.append(get_decyear(str(int(item))));
				dtarrays.append(dtarray);
				decyeararrays.append(decyeararray);
				coords.append([Elong[0]-360, Nlat[0]]);	
				dNs.append(dN);
				dEs.append(dE);
				dUs.append(dU);
				Sns.append(Sn);
				Ses.append(Se);
				Sus.append(Su);
	
	myData=data_collection(name=names,coords=coords, yyyymmdd=yrstrings, dtarray=dtarrays, decyeararray=decyeararrays, dN=dNs, dE=dEs, dU=dUs, Sn=Sns, Se=Ses, Su=Sus);
	return myData;


# ------------------ OUTPUTS -------------------- # 
def outputs(myParams,myData):

	plt.figure();
	[f,axarr]=plt.subplots(1,2,sharex=True,sharey=True)
	label_date="20171031";
	offset=0;
	for i in range(len(myData.dE)):
		offset=16*i;

		[decyear, data_detrend] = detrend_1col_data(myData.decyeararray[i], myData.dE[i]*1000);  # detrend the data in mm. 
		data_offset=[x + offset for x in data_detrend];

		l1 = axarr[0].plot_date(myData.dtarray[i],data_offset,marker='+',markersize=2);
		line_color=l1[0].get_color()
		axarr[0].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,myData.name[i],fontsize=9,color=line_color);
	axarr[0].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20171020", "%Y%m%d"));
	axarr[0].set_ylim([-10,offset+10])
	axarr[0].set_ylabel("East (mm)");
	axarr[0].set_title("Detrended GPS Time Series")
	axarr[0].grid('on')

	for i in range(len(myData.dN)):
		offset=16*i;
		
		[decyear, data_detrend] = detrend_1col_data(myData.decyeararray[i], myData.dN[i]*1000);  # detrend the data in mm. 
		data_offset=[x + offset for x in data_detrend];

		l1 = axarr[1].plot_date(myData.dtarray[i],data_offset,marker='+',markersize=2);
		line_color=l1[0].get_color()
		axarr[1].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,myData.name[i],fontsize=9,color=line_color);
	axarr[1].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20171020", "%Y%m%d"));
	axarr[1].set_ylim([-10,offset+10])
	axarr[1].set_ylabel("North (mm)");
	axarr[1].set_title("Detrended GPS Time Series")
	axarr[1].grid('on')
	plt.savefig(myParams.output_dir+'Collective_TS.jpg')	


	plt.figure();
	for i in range(len(myData.dE)):
		plt.plot(myData.coords[i][0],myData.coords[i][1],'.');
		plt.text(myData.coords[i][0],myData.coords[i][1],myData.name[i]);
	plt.savefig("map.jpg");


	return;



def detrend_1col_data(decyear, data):
	# Detrend the data with a best-fit y = mx + b. 
	coef = np.polyfit(decyear, data, 1);
	data_detrend = [];
	for i in range(len(decyear)):
		data_detrend.append(data[i] - (coef[1]+coef[0]*decyear[i]));
	return [decyear, data_detrend];


def get_decyear(yyyymmdd):
	myday=dt.datetime.strptime(yyyymmdd,"%Y%m%d");
	year=float(yyyymmdd[0:4]);
	jday=float(myday.strftime("%j"));
	decyear=year+jday/365.24;
	return decyear;




