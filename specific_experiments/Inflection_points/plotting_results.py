# Plotting the time of inflection
# Reads the output file from driver.py
# Makes a basemap image. 
# Colorscale is by days since earthquake

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime as dt 
import matplotlib
import matplotlib.cm as cm
from subprocess import call
import gps_ts_functions

def onset_time_map_GMT(name,lon,lat,east,east_change,north, north_change, vert, vert_change, type_of_inflection, earthquake_time):
	# Earthquake time
	eq_dt = dt.datetime.strptime(earthquake_time,"%Y%m%d");
	max_size=3.5;
	outfile_name="Outputs_"+type_of_inflection+"/"+earthquake_time+".txt";
	outfile=open(outfile_name,'w');
	for i in range(len(lon)):
		east_days=(east[i]-eq_dt).days;
		north_days=(north[i]-eq_dt).days;
		vert_days=(vert[i]-eq_dt).days;
		east_size=np.min([abs(east_change[i]),max_size]);
		north_size=np.min([abs(north_change[i]),max_size]);
		vert_size=np.min([abs(vert_change[i]),max_size*1.3]);
		outfile.write("%f %f %f %f %f %f %f %f %f %f %f %s\n" % (lon[i], lat[i], east_change[i], east_size, east_days, north_change[i], north_size, north_days, vert_change[i], vert_size, vert_days, name[i]) );
	outfile.close();
	# output format: lon, lat, east_size, east_days, north_size, north_days, vert_size, vert_days, name

	call(['./timing_map_gps.gmt',outfile_name,'-125', '-118', '36.5', '42.0',type_of_inflection+"/"+earthquake_time+"_norcal",earthquake_time],shell=False);  # northern California
	call(['./timing_map_gps.gmt',outfile_name,'-125', '-110', '32.5', '48.5',type_of_inflection+"/"+earthquake_time+"_WUS",earthquake_time],shell=False);   # WUS
	call(['./timing_map_gps.gmt',outfile_name,'-122', '-115', '32.5', '37.5',type_of_inflection+"/"+earthquake_time+"_socal", earthquake_time],shell=False);   # SoCal
	return;


if __name__=="__main__":

	#  THE MAIN PROGRAM
	# earthquake_time="20140310"  #
	earthquake_time="20161208";
	# type_of_inflection='min_slope';
	type_of_inflection='max_curve';
	infile="Outputs_"+type_of_inflection+"/"+earthquake_time+"_inflections.txt"
	name=[]; lat=[]; lon=[]; east=[]; north=[]; up=[]; east_change=[]; north_change=[]; up_change=[];

	ifile=open(infile,'r');
	for line in ifile:
		temp=line.split();
		name.append(temp[0]);
		lon.append(float(temp[1]));
		lat.append(float(temp[2]));
		east.append(dt.datetime.strptime(temp[3],"%Y-%m-%d"));
		north.append(dt.datetime.strptime(temp[5],"%Y-%m-%d"));
		up.append(dt.datetime.strptime(temp[7],"%Y-%m-%d"));
		east_change.append(float(temp[9]));
		north_change.append(float(temp[10]));
		up_change.append(float(temp[11]));

	# Make GMT maps. 
	onset_time_map_GMT(name, lon, lat, east, east_change, north, north_change, up, up_change, type_of_inflection, earthquake_time);

