# 2/28/2019: The purpose of this plot is to extract the annual amplitudes in GPS time series
# We focus on the vertical compontent first
# And then compare with the equivalent GRACE model that was pre-computed, sitting in a directory somewhere. 

import numpy as np 
import matplotlib.pyplot as plt 
import os
import stations_within_radius
import gps_ts_functions
import gps_input_pipeline
import grace_ts_functions
import offsets
import sys


def configure():
	datasource='pbo';
	refframe='NA';
	grace_dir = "../../GPS_POS_DATA/GRACE_loading_model/"
	coord_box = [-125, -118, 37, 44];
	station_list = stations_within_radius.get_stations_within_box(coord_box, datasource);
	station_list = gps_input_pipeline.remove_blacklist(station_list);
	ampfile = "grace_vs_gps_amps_"+datasource+".txt";
	return [station_list, datasource, refframe, grace_dir, ampfile];

def inputs(station_list, datasource, refframe, grace_dir):
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[]; graceobj_list=[]; 
	for station_name in station_list:
		filename=grace_dir+"scaled_"+station_name+"_PREM_model_ts.txt";
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, datasource, refframe);
		if os.path.isfile(filename):
			print("Opening GRACE for station %s" % station_name);
			grace_obj=grace_ts_functions.input_GRACE_individual_station(filename);
		else:
			print("Error! GRACE not found for %s" % station_name);
			grace_obj=[];
		if myData==[]:
			continue;
		else:
			dataobj_list.append(myData);
			offsetobj_list.append(offset_obj);
			eqobj_list.append(eq_obj);	
			graceobj_list.append(grace_obj);
	return [dataobj_list, offsetobj_list, eqobj_list, graceobj_list];

def compute(dataobj_list, offsetobj_list, eqobj_list, graceobj_list):
	gps_amp=[]; grace_amp=[]; 
	for i in range(len(dataobj_list)):
		print(dataobj_list[i].name);
		newobj=offsets.remove_offsets(dataobj_list[i],offsetobj_list[i]); #  removing earthquakes and offsets
		newobj=offsets.remove_offsets(newobj, eqobj_list[i]);
		if len(newobj.dtarray)<365: # skip the campaign stations
			gps_amp.append([np.nan, np.nan, np.nan]);
			grace_amp.append([np.nan, np.nan, np.nan]);
			continue;
		
		# Get the GPS amplitude
		# Parameter order is slope, cos(wt), sin(wt), cos(2wt), sin(2wt)
		[east_params_gps, north_params_gps, vert_params_gps] = gps_ts_functions.get_linear_annual_semiannual(newobj);	
		east_amp_gps = np.sqrt(east_params_gps[1]*east_params_gps[1]+east_params_gps[2]*east_params_gps[2]);
		north_amp_gps = np.sqrt(north_params_gps[1]*north_params_gps[1]+north_params_gps[2]*north_params_gps[2]);
		vert_amp_gps = np.sqrt(vert_params_gps[1]*vert_params_gps[1]+vert_params_gps[2]*vert_params_gps[2]);

		# Get the GRACE amplitude (skip if we haven't computed the GRACE function yet.)
		if graceobj_list[i]==[]:
			east_amp_grace=np.nan; north_amp_grace=np.nan; vert_amp_grace=np.nan;
		else:
			[east_params_grace, north_params_grace, vert_params_grace] = grace_ts_functions.get_linear_annual_semiannual(graceobj_list[i]);
			east_amp_grace = np.sqrt(east_params_grace[1]*east_params_grace[1]+east_params_grace[2]*east_params_grace[2]);
			north_amp_grace = np.sqrt(north_params_grace[1]*north_params_grace[1]+north_params_grace[2]*north_params_grace[2]);
			vert_amp_grace = np.sqrt(vert_params_grace[1]*vert_params_grace[1]+vert_params_grace[2]*vert_params_grace[2]);
			
		gps_amp.append([east_amp_gps, north_amp_gps, vert_amp_gps]);
		grace_amp.append([east_amp_grace, north_amp_grace, vert_amp_grace]);
	return [gps_amp, grace_amp];

def outputs(dataobj_list, gps_amp, grace_amp, ampfile):
	ofile=open(ampfile,'w');
	for i in range(len(dataobj_list)):
		ofile.write("%f %f %f %f %f %f %f %f %s\n" % (dataobj_list[i].coords[0], dataobj_list[i].coords[1], gps_amp[i][0], gps_amp[i][1], gps_amp[i][2], grace_amp[i][0], grace_amp[i][1], grace_amp[i][2], dataobj_list[i].name) );
	ofile.close();
	return;

def get_percent_close(gps_amp, grace_amp, close_limit):
	close_count = 0;
	total = 0;
	for i in range(len(gps_amp)):
		if abs(gps_amp[i] - grace_amp[i])<=close_limit:
			close_count=close_count+1;
		if abs(gps_amp[i] - grace_amp[i])>=0:  # for filtering out the nans
			total = total+1;
	percent_close = 100* (close_count/total);
	return percent_close;

def remake_plots():
	# Make a pretty plot 
	unrfile="grace_vs_gps_amps_unr.txt"
	pbofile="grace_vs_gps_amps_pbo.txt"
	[lon_unr, lat_unr, east_gps_unr, gps_amp_unr, east_grace_unr, grace_amp_unr]=np.loadtxt(unrfile,unpack=True,usecols=[0,1,2,4,5,7]);
	[lon_pbo, lat_pbo, east_gps_pbo, gps_amp_pbo, east_grace_pbo, grace_amp_pbo]=np.loadtxt(pbofile,unpack=True,usecols=[0,1,2,4,5,7]);

	# Defining a metric of how many stations have GPS and GRACE amplitudes that are pretty close. 
	close_limit = 1;  # mm
	unr_percent = get_percent_close(gps_amp_unr, grace_amp_unr, close_limit);
	pbo_percent = get_percent_close(gps_amp_pbo, grace_amp_pbo, close_limit);

	# The plotting guts. 
	fig = plt.figure();
	h1 = plt.plot(gps_amp_unr, grace_amp_unr,'.',markersize=4,label='UNR: %d%% close to matching' % (unr_percent) );
	h2 = plt.plot(gps_amp_pbo, grace_amp_pbo,'.',markersize=4,label='PBO: %d%% close to matching' % (pbo_percent) );
	plt.xlabel('GPS Seasonal Amplitude (mm)');
	plt.ylabel('GRACE Seasonal Amplitude (mm)');
	mmax=10;
	plt.plot([0,mmax],[0,mmax],'--k');
	plt.plot([0,mmax],[0-close_limit,mmax-close_limit],'--',color='gray');
	plt.plot([0,mmax],[0+close_limit,mmax+close_limit],'--',color='gray');
	plt.xlim([0,mmax]);
	plt.ylim([0,mmax]);
	plt.legend()
	plt.savefig('vert_amp_vs_amp_both.eps');


	# Defining a metric of how many stations have GPS and GRACE amplitudes that are pretty close. 
	close_limit = 1;  # mm
	unr_percent = get_percent_close(east_gps_unr, east_grace_unr, close_limit);
	pbo_percent = get_percent_close(east_gps_pbo, east_grace_pbo, close_limit);	

	fig = plt.figure();
	h1 = plt.plot(east_gps_unr, east_grace_unr,'.',markersize=4,label='UNR: %d%% close to matching' % (unr_percent) );
	h2 = plt.plot(east_gps_pbo, east_grace_pbo,'.',markersize=4,label='PBO: %d%% close to matching' % (pbo_percent) );
	plt.xlabel('GPS Seasonal Amplitude (mm)');
	plt.ylabel('GRACE Seasonal Amplitude (mm)');
	mmax=3;
	plt.plot([0,mmax],[0,mmax],'--k');
	plt.plot([0,mmax],[0-close_limit,mmax-close_limit],'--',color='gray');
	plt.plot([0,mmax],[0+close_limit,mmax+close_limit],'--',color='gray');
	plt.xlim([0,mmax]);
	plt.ylim([0,mmax]);
	plt.legend()
	plt.savefig('east_amp_vs_amp_both.eps');	

	return;


if __name__=="__main__":
	# [station_list, datasource, refframe, grace_dir, ampfile] = configure();
	# [dataobj_list, offsetobj_list, eqobj_list, graceobj_list] = inputs(station_list, datasource, refframe, grace_dir);
	# [gps_amp, grace_amp] = compute(dataobj_list, offsetobj_list, eqobj_list, graceobj_list);
	# outputs(dataobj_list, gps_amp, grace_amp, ampfile);

	remake_plots();

