# Python script to download a bunch of UNR Data
# Data source: http://geodesy.unr.edu/PlugNPlayPortal.php
# ref_frame is either ISG08 or NA12.  (IGS14 as of November 2019)
# July 2018
# For most of California, it's about 1200 stations and about 500Mb. It takes ~10 minutes. 
# For the western US, it's about 2600 stations and ~20 minutes. 
# One hour if the data transfer is slow. 

import numpy as np 
import sys
from subprocess import call

def configure():
	coordfile = "UNR_coords_dec2019.txt";  # the data holding file (for getting our region of interest)
	ref_frame1 = "NA";
	ref_frame2 = "IGS08";
	ref_frame3 = "IGS14";
	latlon_box = [-125.0, -110, 32.0, 49.0];  # A LARGE BOX INCLUDING ALL OF WUS
	return coordfile, ref_frame1, ref_frame2, ref_frame3, latlon_box;

def get_stations(coordfile, latlon_box):
	station_names=[];
	data = np.genfromtxt(coordfile,skip_header=2,usecols=(0,1,2), dtype=None);
	for item in data:
		lat=item[1];
		lon=item[2];
		if lon>180:
			lon=lon-360;
		if lon>=latlon_box[0] and lon<=latlon_box[1]:
			if lat>=latlon_box[2] and lat<=latlon_box[3]:
				onestation = item[0].decode('UTF-8');
				station_names.append(onestation);
	print("Returning %d stations for downloading. " % len(station_names));
	return station_names;

def download_stations(stations, ref_frame):
	if ref_frame=="NA":
		ref_frame_local="NA";
	else:
		ref_frame_local=ref_frame;
	for mystation in stations:
		print("Downloading station %s " % mystation);
		if ref_frame=="IGS14":
			call('wget http://geodesy.unr.edu/gps_timeseries/tenv3/'+ref_frame+'/'+mystation+'.tenv3 -O '+mystation+'.'+ref_frame+'.tenv3',shell=True);
		else:
			call('wget http://geodesy.unr.edu/gps_timeseries/tenv3/'+'plates'+'/'+ref_frame+'/'+mystation+'.'+ref_frame+'.tenv3 -O '+mystation+'.'+ref_frame_local+'.tenv3',shell=True);
	return;


if __name__=="__main__":
	coordfile, ref_frame1, ref_frame2, ref_frame3, latlon_box = configure();
	stations = get_stations(coordfile, latlon_box);
	download_stations(stations, ref_frame1);
	# download_stations(stations, ref_frame1);
	# download_stations(stations, ref_frame2);
