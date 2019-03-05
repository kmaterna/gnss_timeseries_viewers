# Read in the 2005 GPS stations
# Calculate the total afterslip displacements. 
# One year is good enough. 


import numpy as np
import matplotlib.pyplot as plt 
import collections
import datetime as dt 
import stations_within_radius
import gps_ts_functions
import gps_seasonal_removals
import offsets
import gps_input_pipeline


def configure():
	box = [-125, -120, 38, 42.1];
	stations = stations_within_radius.get_stations_within_box(box);
	proc_center='pbo';
	early_date = dt.datetime.strptime("20050616","%Y%m%d");
	later_date = dt.datetime.strptime("20060615","%Y%m%d");
	filename="2005_gps_disps.txt"
	return [stations, proc_center, early_date, later_date, filename]; 

# INPUTS
def inputs(station_names, proc_center):  # Returns a list of objects for time series data, offsets, and earthquakes
	dataobj_list=[]; offsetobj_list=[]; eqobj_list=[];
	for station_name in station_names:
		[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station_name, proc_center, "NA");
		if myData.dtarray[0]<=dt.datetime.strptime("20050615","%Y%m%d"):
		# kicking out the stations that end early or start late. 
			dataobj_list.append(myData);
			offsetobj_list.append(offset_obj);
			eqobj_list.append(eq_obj);	
	return [dataobj_list, offsetobj_list, eqobj_list];

# COMPUTE
def compute(dataobj_list, offsetobj_list, eqobj_list, early_date, later_date):
	east_afterslip = []; north_afterslip = []; 
	for i in range(len(dataobj_list)):
		# Remove the steps from earthquakes and seasonals. 
		newobj=offsets.remove_offsets(dataobj_list[i], offsetobj_list[i]);
		newobj=offsets.remove_offsets(newobj,eqobj_list[i]);
		newobj=gps_seasonal_removals.make_detrended_ts(newobj, 1, 'lssq');
		# compute afterslip displacement
		east_disp = get_displacements(newobj.dtarray, newobj.dE, early_date, later_date);
		north_disp = get_displacements(newobj.dtarray, newobj.dN, early_date, later_date);
		east_afterslip.append(east_disp);
		north_afterslip.append(north_disp);
	return [east_afterslip, north_afterslip];


def get_displacements(t, x, early_time, late_time):
	first_measurement=np.nan;
	last_measurement=np.nan;
	for i in range(len(t)):
		if t[i]==early_time:
			first_measurement = np.mean(x[i:i+2]);
		elif t[i]==late_time:
			last_measurement = np.mean(x[i:i+5]);  # averaging helps avoid noise 
		else:
			continue;
	return last_measurement - first_measurement;

def write_outputs(dataobj_list, east_afterslip, north_afterslip, filename):
	ofile=open(filename,'w');
	for i in range(len(dataobj_list)):
		ofile.write("%f %f %f %f 0 0 0 %s\n" % (dataobj_list[i].coords[0], dataobj_list[i].coords[1], east_afterslip[i], north_afterslip[i], dataobj_list[i].name) );
	ofile.close();
	return;

if __name__=="__main__":
	[stations, proc_center, early_date, later_date, filename] = configure();
	[dataobj_list, offsetobj_list, eqobj_list] = inputs(stations, proc_center);
	[east_afterslip, north_afterslip] = compute(dataobj_list, offsetobj_list, eqobj_list, early_date, later_date);
	write_outputs(dataobj_list, east_afterslip, north_afterslip, filename);

