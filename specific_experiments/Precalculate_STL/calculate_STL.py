# The purpose of this script is to perform STL on each station in the list. 
# It will save the data into a text file
# This process will save time down the line. 
# 30 may not always be the optimal parameter value (especially for long time series), 
# But I'm leaving it hard-coded for now. 

import numpy as np 
import collections
import datetime as dt
import gps_input_pipeline
import offsets
import gps_seasonal_removals
import gps_ts_functions

Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

def main_function():
	[infile,outdir, outliers_def]=configure();
	[DataAll, OffsetAll, EQAll]=inputs(infile);
	compute(DataAll,OffsetAll, EQAll, outliers_def, outdir);
	return;

def configure():
	infile="CA2.txt";
	outdir="../../GPS_POS_DATA/STL_models/";
	outliers_def = 15.0;  # mm away from average. 
	return [infile, outdir, outliers_def];

def inputs(infile):
	DataAll=[]; OffsetAll=[]; EQAll=[];
	ifile=open(infile,'r');
	for line in ifile:
		temp=line.split();
		station_name=temp[0];
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name,'pbo');
		DataAll.append(myData);
		OffsetAll.append(offset_obj);
		EQAll.append(eq_obj);
	ifile.close();
	return [DataAll, OffsetAll, EQAll];


# -------------- COMPUTE ------------ # 
def compute(myData, offset_obj, eq_obj, outliers_def, outdir):
	for i in range(len(myData)):
		newData=myData[i]; 
		newData=offsets.remove_antenna_offsets(newData, offset_obj[i]);
		newData=gps_ts_functions.remove_outliers(newData, outliers_def);
		newData=offsets.remove_earthquakes(newData, eq_obj[i]);
		stl_filt=gps_seasonal_removals.make_detrended_ts(newData, 1, 'stl');
		outputs(stl_filt, outdir);
	return;

# ------------ OUTPUTS ------------- # 
def outputs(Data, outdir):
	ofile=open(outdir+Data.name+"_STL_30.txt",'w');
	for i in range(len(Data.dtarray)):
		timestamp=dt.datetime.strftime(Data.dtarray[i],"%Y%m%d");
		ofile.write("%s %f %f %f %f %f %f\n" % (timestamp, Data.dE[i], Data.dN[i], Data.dU[i], Data.Se[i], Data.Sn[i], Data.Su[i]) );
	ofile.close();
	return;


if __name__=="__main__":
	main_function();