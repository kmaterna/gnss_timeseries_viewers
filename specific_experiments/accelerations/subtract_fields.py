# Subtract two datasets

import numpy as np 
import collections
import subprocess

Velfield=collections.namedtuple("Velfield",['lon','lat','east','north','vert','name']);

def configure():
	year="2016";
	network="UNR";
	file1="Whole_WUS_lssq_"+network+"_nam/MTJ_"+year+"_lssq.txt";
	file2="Whole_WUS_lssq_"+network+"_itrf/MTJ_"+year+"_lssq.txt";
	outfile="subtracted.txt";
	return [file1, file2, outfile];

def input_field(input_file):
	lons=[]; lats=[]; east=[]; north=[]; vert=[]; name=[];
	ifile=open(input_file,'r');
	for line in ifile:
		temp=line.split();
		lons.append(float(temp[0]));
		lats.append(float(temp[1]));
		east.append(float(temp[2]));
		north.append(float(temp[3]));
		vert.append(float(temp[5]));
		name.append(temp[8]);
	ifile.close();
	myvels=Velfield(lon=lons, lat=lats, east=east, north=north, vert=vert, name=name);
	return myvels;

def subtract(vel1, vel2):
	name=[]; lon=[]; lat=[]; east=[]; north=[]; vert=[];
	for i in range(len(vel1.lon)):
		if vel1.name[i] in vel2.name:
			idx=vel2.name.index(vel1.name[i]);
			name.append(vel1.name[i]);
			lon.append(vel1.lon[i]);
			lat.append(vel1.lat[i]);
			east.append(vel1.east[i]-vel2.east[idx]);
			north.append(vel1.north[i]-vel2.north[idx]);
			vert.append(vel1.vert[i]-vel2.vert[idx]);
	Diffs=Velfield(lon=lon, lat=lat, east=east, north=north, vert=vert, name=name);
	return Diffs;

def outputs(Diffs, outfile):
	ofile=open(outfile,'w');
	for i in range(len(Diffs.name)):
		ofile.write("%f %f %f %f 0 %f 0 0 %s\n" % (Diffs.lon[i], Diffs.lat[i], Diffs.east[i], Diffs.north[i], Diffs.vert[i], Diffs.name[i]) );
	ofile.close();
	subprocess.call('./accel_map_gps.gmt '+outfile+' -125.2 -121.0 38.6 43.0 MTJ_differences',shell=True);
	subprocess.call('./accel_map_gps.gmt '+outfile+' -124.6 -118.4 35.5 42.2 NorCal_differences',shell=True);
	subprocess.call('./accel_map_gps.gmt '+outfile+' -121.8 -115.0 32.2 37.6 SoCal_differences',shell=True);
	subprocess.call('./accel_map_gps.gmt '+outfile+' -125.6 -110.0 32.5 48.5 WUS_differences',shell=True);
	subprocess.call('./accel_map_gps.gmt '+outfile+' -123.5 -119.0 35.6 40.0 SF_differences',shell=True);
	subprocess.call('./accel_map_gps.gmt '+outfile+' -124.6 -120.4 41.2 46.2 Oregon_differences',shell=True);
	return;

if __name__=="__main__":
	[file1, file2, outfile]=configure();
	vel1=input_field(file1);
	vel2=input_field(file2);
	Diffs=subtract(vel1,vel2);
	outputs(Diffs, outfile);

