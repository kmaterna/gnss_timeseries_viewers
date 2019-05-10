# Python script
# April 21, 2019
# This script takes a file of stations (GRACE computation compatible)
# And downloads the loading products for each station from the German loading website
# http://rz-vm115.gfz-potsdam.de:8080/repository/entry/show?entryid=362f8705-4b87-48d1-9d86-2cfd1a2b6ac9
# I use the Center of Figure HYDL products. 

import sys
import subprocess
import datetime as dt 
import gps_ts_functions

def configure():
	input_file="norcal_subset.txt";
	output_dir="../../../GPS_POS_DATA/PBO_Hydro/LSDM/";
	return [input_file, output_dir];

def get_stations(input_file, output_dir):
	ifile=open(input_file,'r');
	for line in ifile:
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

		# The actual calling of the script
		if station_name=="CME6":
			continue;
		print("Getting LSDM Hydro for station %s " % station_name);
		print("extractlatlon_bilinintp_remote HYDL CF "+dt.datetime.strftime(start_date,"%Y-%m-%d")+" "+dt.datetime.strftime(end_date,"%Y-%m-%d")+" "+lat+" "+lon+" -o "+output_dir+station_name+"_LSDM_hydro.txt");
		subprocess.call(["./extractlatlon_bilinintp_remote HYDL CF "+dt.datetime.strftime(start_date,"%Y-%m-%d")+" "+dt.datetime.strftime(end_date,"%Y-%m-%d")+" "+lat+" "+lon+" -o "+output_dir+station_name+"_LSDM_hydro.txt"],shell=True);
	ifile.close();
	return;


if __name__=="__main__":
	[input_file, output_dir]=configure();
	get_stations(input_file, output_dir);