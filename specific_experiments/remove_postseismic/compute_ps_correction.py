# Remove postseismic from Salton Sea GNSS stations using Hines et al., JGR, 2016 model. 
# -------- STEPS -------- # 
# Read trever's file into a list of TS objects
# Interpolate over time, simple 1D
# Interpolate over space too, using NP 2D interp (a triangulation method) for needed stations
# Write a set of .pos files for each required station
# This is meant to be done only once in a while, since we can use the modeled TS afterwards for correction. 


import numpy as np 
import h5py
import matplotlib.pyplot as plt 
import datetime as dt
import sys
from scipy import interpolate
import gps_ts_functions
import gps_io_functions
import stations_within_radius

# --------------- HELPER FUNCTIONS --------------- # 

def interpolate_over_time(dates, east, north, vert):
	# This function will take a set of dates and associated smooth functions for east, north, and up displacements
	# and interpolate the smooth functions into daily predictions.
	dates_intended = gps_ts_functions.get_daily_dtarray(dates[0], dates[-1]);
	floattimes_actual = gps_ts_functions.get_float_times(dates);
	floattimes_intended = gps_ts_functions.get_float_times(dates_intended);
	f = interpolate.interp1d(floattimes_actual, east);
	y_east = f(floattimes_intended);
	f = interpolate.interp1d(floattimes_actual, north);
	y_north = f(floattimes_intended);
	f = interpolate.interp1d(floattimes_actual, vert);
	y_vert = f(floattimes_intended);	
	return dates_intended, y_east, y_north, y_vert;


# --------------- THE PROGRAM --------------- # 

def configure():
	intended_box = [-116, -115, 32.5, 33.5];
	intended_stations1 = stations_within_radius.get_stations_within_box(intended_box, network='unr');
	intended_stations2 = stations_within_radius.get_stations_within_box(intended_box, network='pbo');
	intended_stations = list(set(intended_stations1+intended_stations2));
	print("Returning %d intended stations " % (len(intended_stations)) );
	model_file = "../../GPS_POS_DATA/Remove_postseismic/Hines/results.h5"
	starttime = dt.datetime.strptime("2010-04-04","%Y-%m-%d");  # earthquake time
	endtime = dt.datetime.strptime("2020-04-04","%Y-%m-%d");  # the end of the time series we're writing
	return intended_stations, model_file, starttime, endtime;

# --------------- INPUTS --------------- # 

def read_hines_to_tsObjList(filename, starttime):
	print("Reading %s into TS objects" % (filename) )
	f = h5py.File(filename,'r');
	print(f.keys())

	tsObjList = [];

	model = f['predicted_displacement'];  # the full viscoelastic + afterslip model
	num_objects = len(model['name']);
	for i in range(num_objects):
		obj_name = model['name'][i].decode('utf-8');
		# print(obj_name);
		time = model['time'];
		coords = model['position'][i];
		mystation = model['mean'][:,i];
		east_pred = [x[0]*1000 for x in mystation];
		north_pred = [x[1]*1000 for x in mystation];
		vert_pred = [x[2]*1000 for x in mystation];
		# tsobjects have displacement in millimeters

		# Make a dtarray that matches the time array, and interpolate the values through that array
		dtarray = gps_ts_functions.yrnum2datetime(time, starttime);
		dates_intended, y_east, y_north, y_vert = interpolate_over_time(dtarray, east_pred, north_pred, vert_pred);

		myts = gps_ts_functions.Timeseries(name=obj_name, coords=coords, dtarray=dates_intended, 
			dN=y_north, dE=y_east, dU=y_vert, 
			Sn=np.zeros(np.shape(y_north)), Se=np.zeros(np.shape(y_east)), Su=np.zeros(np.shape(y_vert)), EQtimes=[]);
		tsObjList.append(myts);

	return tsObjList;

def how_many_new_stations(intended_stations, tsObjList):
	given_stations = [tsObjList[i].name for i in range(len(tsObjList))];
	new_stations = []; existing_stations=[]; 
	for i in range(len(intended_stations)):
		if intended_stations[i] not in given_stations:
			new_stations.append(intended_stations[i]);
		else:
			existing_stations.append(intended_stations[i]);
	[new_lon, new_lat] = gps_io_functions.get_coordinates_for_stations(new_stations);
	print("We have modeled timeseries for %d desired stations " % (len(existing_stations)) )
	print("Need to compute for %d new stations " % len(new_stations));
	return new_stations, new_lon, new_lat;


# --------------- COMPUTE --------------- # 

def compute_for_new_stations(new_stations, new_lon, new_lat, tsObjList):
	# Perform a 2D interpolation for each day. 
	print("Interpolating into new stations"); 
	newTsObjList=[];
	new_dE=[]; new_dN = []; new_dU = [];  # lists of lists
	for item in new_stations:
		new_dE.append([]); new_dN.append([]); new_dU.append([]);

	unpack_lon = [tsObjList[x].coords[0] for x in range(len(tsObjList))];
	unpack_lat = [tsObjList[x].coords[1] for x in range(len(tsObjList))];
	given_coords = list(zip(unpack_lon, unpack_lat));

	# For each day in the time arrays... 
	for i in range(len(tsObjList[0].dtarray)):
		unpack_time = [tsObjList[x].dtarray[i] for x in range(len(tsObjList))];
		unpack_east = [tsObjList[x].dE[i] for x in range(len(tsObjList))];
		unpack_north = [tsObjList[x].dN[i] for x in range(len(tsObjList))];
		unpack_vert = [tsObjList[x].dU[i] for x in range(len(tsObjList))];

		# Scipy returns a function that you can use on a new set of x,y pairs. 
		f_east = interpolate.LinearNDInterpolator(given_coords, unpack_east, fill_value=0);
		f_north = interpolate.LinearNDInterpolator(given_coords, unpack_north, fill_value=0);
		f_vert = interpolate.LinearNDInterpolator(given_coords, unpack_vert, fill_value=0);

		for j in range(len(new_stations)):
			# Compute the interpolation for each required station
			new_east=f_east(new_lon[j], new_lat[j]);  # only want to give the functions one point at a time. 
			new_north=f_north(new_lon[j], new_lat[j]);
			new_vert =f_vert(new_lon[j], new_lat[j]);

			# Pack the results back up. 
			new_dE[j].append(new_east);
			new_dN[j].append(new_north);
			new_dU[j].append(new_vert);

	for j in range(len(new_stations)):
		myts = gps_ts_functions.Timeseries(name=new_stations[j], coords=[new_lon[j], new_lat[j]], dtarray=tsObjList[0].dtarray, 
			dN=new_dN[j], dE=new_dE[j], dU=new_dU[j], 
			Sn=np.zeros(np.shape(new_dN[j])), Se=np.zeros(np.shape(new_dE[j])), Su=np.zeros(np.shape(new_dU[j])), EQtimes=[]);
		newTsObjList.append(myts);

	return newTsObjList;

def lengthen_timeseries(tsObjList, starttime, endtime):
	# Using a log model to fit the time series and lengthen it. 
	newTsObjList=[];
	for i in range(len(tsObjList)):

		Data0 = tsObjList[i];
		[eparams, nparams, uparams] = gps_ts_functions.get_logfunction(Data0, starttime);

		longer_ts = gps_ts_functions.get_daily_dtarray(starttime, endtime);
		longer_float_time = gps_ts_functions.get_relative_times(longer_ts, starttime);

		e_model = gps_ts_functions.construct_log_function(longer_float_time, eparams);
		n_model = gps_ts_functions.construct_log_function(longer_float_time, nparams);
		u_model = gps_ts_functions.construct_log_function(longer_float_time, uparams);

		# This is a little tricky.  Here we are appending the extrapolated log fit
		# to the existing Hines' model. 
		new_e = list(Data0.dE); new_n = list(Data0.dN); new_u = list(Data0.dU);
		for i in range(len(longer_ts)):
			if longer_ts[i]>Data0.dtarray[-1]:
				new_e.append(e_model[i]);
				new_n.append(n_model[i]);
				new_u.append(u_model[i]);

		myts = gps_ts_functions.Timeseries(name=Data0.name, coords=Data0.coords, dtarray=longer_ts, 
			dE=new_e, dN=new_n, dU=new_u, 
			Sn=np.zeros(np.shape(new_n)), Se=np.zeros(np.shape(new_e)), Su=np.zeros(np.shape(new_u)), EQtimes=[]);

		newTsObjList.append(myts);
		
	return newTsObjList;



# --------------- OUTPUTS --------------- # 

def write_model_ts(tsObjList):
	for x in tsObjList:
		filename="../../GPS_POS_DATA/Remove_postseismic/Hines/Stations/"+x.name+"_psmodel.pos";
		gps_io_functions.write_pbo_pos_file(x, filename);
	return;

def interpolation_map_figure(tsObjList, newTsObjList):

	# Unpacking
	new_lon=[x.coords[0] for x in newTsObjList];
	new_lat=[x.coords[1] for x in newTsObjList];
	xlon=[x.coords[0] for x in tsObjList];
	ylat=[x.coords[1] for x in tsObjList];
	dE_modeled= [x.dE[-1]-x.dE[0] for x in tsObjList];
	dE_interp = [x.dE[-1]-x.dE[0] for x in newTsObjList];
	dN_modeled= [x.dN[-1]-x.dN[0] for x in tsObjList];
	dN_interp = [x.dN[-1]-x.dN[0] for x in newTsObjList];
	dU_modeled= [x.dU[-1]-x.dU[0] for x in tsObjList];
	dU_interp = [x.dU[-1]-x.dU[0] for x in newTsObjList];
	[ca_lon, ca_lat] = np.loadtxt('california_bdr',unpack=True);
	[az_lon, az_lat] = np.loadtxt('arizona_bdr',unpack=True);

	plt.figure(dpi=300);
	plt.plot(ca_lon,ca_lat,'k');
	plt.plot(az_lon,az_lat,'k');
	plt.scatter(xlon, ylat, s=12, c=dE_modeled, marker='s', cmap='jet', vmin=-10, vmax=10);
	plt.scatter(new_lon, new_lat, s=23, c=dE_interp, cmap='jet', vmin=-10, vmax=10);
	plt.xlim([-117, -114]);
	plt.ylim([32.4, 34.5]);
	plt.title('Interpolated East (circles) with modeled (squares)');
	plt.ylabel('Latitude',fontsize=15);
	plt.xlabel('Longitude',fontsize=15);
	cbar = plt.colorbar();
	cbar.ax.set_ylabel('5-year Displacement (mm)');
	plt.savefig("interpolated_east_displacements.png");

	plt.figure(dpi=300);
	plt.plot(ca_lon,ca_lat,'k');
	plt.plot(az_lon,az_lat,'k');
	plt.scatter(xlon, ylat, s=12, c=dN_modeled, marker='s', cmap='jet', vmin=-40, vmax=10);
	plt.scatter(new_lon, new_lat, s=23, c=dN_interp, cmap='jet', vmin=-40, vmax=10);
	plt.xlim([-117, -114]);
	plt.ylim([32.4, 34.5]);
	plt.title('Interpolated North (circles) with modeled (squares)');
	plt.ylabel('Latitude',fontsize=15);
	plt.xlabel('Longitude',fontsize=15);	
	cbar = plt.colorbar();
	cbar.ax.set_ylabel('5-year Displacement (mm)');
	plt.savefig("interpolated_north_displacements.png");	

	plt.figure(dpi=300);
	plt.plot(ca_lon,ca_lat,'k');
	plt.plot(az_lon,az_lat,'k');
	plt.scatter(xlon, ylat, s=12, c=dU_modeled, marker='s', cmap='jet', vmin=-15, vmax=15);
	plt.scatter(new_lon, new_lat, s=23, c=dU_interp, cmap='jet', vmin=-15, vmax=15);
	plt.xlim([-117, -114]);
	plt.ylim([32.4, 34.5]);
	plt.title('Interpolated Up (circles) with modeled (squares)');
	plt.ylabel('Latitude',fontsize=15);
	plt.xlabel('Longitude',fontsize=15);	
	cbar = plt.colorbar();
	cbar.ax.set_ylabel('5-year Displacement (mm)');
	plt.savefig("interpolated_up_displacements.png");

	return;

def timeseries_figure(ObjList):
	# A 3-panel plot with all the modeled time series
	# It might be good for the supplement. 
	fig, axarr = plt.subplots(3,1, figsize=(6,12));

	for x in ObjList:
		axarr[0].plot(x.dtarray, x.dE);
		axarr[1].plot(x.dtarray, x.dN);
		axarr[2].plot(x.dtarray, x.dU);

	# Formatting
	axarr[0].set_ylabel('East (mm)');
	axarr[0].set_title('Modeled Time Series from Hines et al., 2016');		
	axarr[1].set_ylabel('North (mm)');
	axarr[2].set_ylabel('Vertical (mm)');
	axarr[2].set_xlabel('Time');

	fig.savefig("Hines_modeled_timeseries.png");

	return;

def view_log_model(tsObjList, modeled_tsObjList):
	# This is a separate function to plot the time series of stations
	# Plus their best-fitting log function (to see how well they match)
	# Might be good in the supplement. 
	plt.figure();
	for i in range(len(tsObjList)):
		Data0 = tsObjList[i];
		Data0_longer = modeled_tsObjList[i];
		plt.plot(Data0.dtarray, Data0.dN,'.');
		plt.plot(Data0_longer.dtarray, Data0_longer.dN, '--k');
	plt.title('Hines Model and Logarithmic Fit: North');
	plt.savefig("log_functions_north.png");

	plt.figure();
	for i in range(len(tsObjList)):
		Data0 = tsObjList[i];
		Data0_longer = modeled_tsObjList[i];
		plt.plot(Data0.dtarray, Data0.dE,'.');
		plt.plot(Data0_longer.dtarray, Data0_longer.dE, '--k');
	plt.title('Hines Model and Logarithmic Fit: East');
	plt.savefig("log_functions_east.png");	

	plt.figure();
	for i in range(len(tsObjList)):
		Data0 = tsObjList[i];
		Data0_longer = modeled_tsObjList[i];
		plt.plot(Data0.dtarray, Data0.dU,'.');
		plt.plot(Data0_longer.dtarray, Data0_longer.dU, '--k');
	plt.title('Hines Model and Logarithmic Fit: Up');
	plt.savefig("log_functions_up.png");		
	return;


if __name__=="__main__":
	# Where are the new Salton Sea stations? Configure and Input
	intended_stations, model_file, starttime, endtime = configure();
	tsObjList = read_hines_to_tsObjList(model_file, starttime);

	# # Outputs directly from the Hines model. 
	# # Useful for supplementary figures. 
	# view_log_model(tsObjList);
	# timeseries_figure(tsObjList);

	# Can we interpolate into new stations? 
	new_stations, new_lon, new_lat = how_many_new_stations(intended_stations, tsObjList);
	newTsObjList = compute_for_new_stations(new_stations, new_lon, new_lat, tsObjList);
	# interpolation_map_figure(tsObjList, newTsObjList);
	
	# Can we lengthen the time series through a model? 
	modeled_tsObjList = lengthen_timeseries(tsObjList, starttime, endtime);
	modeled_newTsObjList = lengthen_timeseries(newTsObjList, starttime, endtime);

	# OUTPUTS
	write_model_ts(modeled_tsObjList);
	write_model_ts(modeled_newTsObjList);
	timeseries_figure(tsObjList);  # useful for supplement
	view_log_model(tsObjList, modeled_tsObjList);  # useful for supplement



