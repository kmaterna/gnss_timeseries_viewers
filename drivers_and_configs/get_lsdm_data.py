# Python script
# April 21, 2019
# This script takes a file of stations (GRACE computation compatible)
# And downloads the loading products for each station from the German loading website
# http://rz-vm115.gfz-potsdam.de:8080/repository/entry/show?entryid=362f8705-4b87-48d1-9d86-2cfd1a2b6ac9
# I use the Center of Figure HYDL products. 
# Put your desired stations into a text file with row format [name longitude latitude startdate enddate]

import sys, os
import subprocess
import datetime as dt 
import gps_ts_functions

def configure():
	input_file="CA_OR.txt";
	output_dir="../../../GPS_POS_DATA/PBO_Hydro/LSDM/";  # or wherever you want to put your LSDM loads
	product="HYDL";  # OPTIONS: HYDL, NTAL, NTOL+NTAL, etc. 
	return [input_file, product, output_dir];

def get_stations(input_file, product, output_dir):
	ifile=open(input_file,'r');
	ofile=open("new_stations.txt",'w');
	ofile.close();
	for line in ifile:
		ofile=open("new_stations.txt",'a');
		temp=line.split();
		if len(line.split())==0:
			continue;
		station_name=temp[0];
		lon=temp[1];
		lat=temp[2];
		dec_startdate=float(temp[3]);
		dec_enddate=float(temp[4]);
		start_date=gps_ts_functions.float_to_dt(dec_startdate);
		end_date=gps_ts_functions.float_to_dt(dec_enddate);

		if os.path.isfile(output_dir+station_name+"_LSDM_hydro.txt.txt"): # if we already have the file
			print("Skipping "+station_name+ " because we already have the file");
			continue;

		print("Getting LSDM Hydro for station %s " % station_name);
		print("extractlatlon_bilinintp_remote "+product+" CF "+dt.datetime.strftime(start_date,"%Y-%m-%d")+" "+dt.datetime.strftime(end_date,"%Y-%m-%d")+" "+lat+" "+lon+" -o "+output_dir+station_name+"_LSDM_hydro.txt");
		subprocess.call(["./extractlatlon_bilinintp_remote "+product+" CF "+dt.datetime.strftime(start_date,"%Y-%m-%d")+" "+dt.datetime.strftime(end_date,"%Y-%m-%d")+" "+lat+" "+lon+" -o "+output_dir+station_name+"_LSDM_hydro.txt"],shell=True);
		ofile.write(line+"\n");
		ofile.close();
	ifile.close();
	return;


if __name__=="__main__":
	[input_file, product, output_dir]=configure();
	get_stations(input_file, product, output_dir);