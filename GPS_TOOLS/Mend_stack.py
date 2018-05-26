# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and eventually trend from the data. 
# Step 4: Plot in order of increasing latitude, colored by how close they are to the earthquake. 

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import matplotlib.pyplot as plt 
import collections
import datetime as dt 
from scipy import signal
import gps_io_functions
import gps_ts_functions
import stations_within_radius


def driver():
	[stations, distances, filenames, EQtime, earthquakes_dir] = configure();
	dataobj_list = inputs(stations, filenames);
	[sorted_objects, noeq_objects, sorted_distances] = compute(dataobj_list, distances, earthquakes_dir);
	
	output_full_ts(sorted_objects, sorted_distances, EQtime, "detrended");
	output_full_ts(noeq_objects, sorted_distances, EQtime, "noeq");
	return;


def configure():
	EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
	EQtime  = dt.datetime.strptime("20140310", "%Y%m%d");
	earthquakes_dir = earthquakes_dir="../GPS_POS_DATA/Event_Files/";
	radius=120;  # km. 
	stations, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius);
	filenames=[];
	for station in stations:
		filenames.append("../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos");
	return [stations, distances, filenames, EQtime, earthquakes_dir];

def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;


def compute(dataobj_list, distances, earthquakes_dir):

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

	# Detrended objects
	detrended_objects=[];
	for i in range(len(dataobj_list)):
		newobj=gps_ts_functions.detrend_data(sorted_objects[i]);
		detrended_objects.append(newobj);


	# No earthquakes objects
	# SOMETHING ELSE HAPPENS
	noeq_objects = [];
	for i in range(len(dataobj_list)):
		newobj=gps_ts_functions.remove_earthquakes(sorted_objects[i],earthquakes_dir);
		newobj=gps_ts_functions.detrend_data(newobj);
		noeq_objects.append(newobj);

	return [detrended_objects, noeq_objects, sorted_distances];


def output_full_ts(dataobj_list, distances, EQtime, filename):

	plt.figure();
	[f,axarr]=plt.subplots(1,2,sharex=True,sharey=True)
	label_date="20171031";
	offset=0;
	spacing=10;
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		edata=dataobj_list[i].dE;
		edata=[x + offset for x in edata];
		l1 = axarr[0].plot_date(dataobj_list[i].dtarray,edata,marker='+',markersize=2);
		line_color=l1[0].get_color()
		axarr[0].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,dataobj_list[i].name,fontsize=9,color=line_color);
	axarr[0].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20171020", "%Y%m%d"));
	axarr[0].set_ylim([-10,offset+10])
	bottom,top=axarr[0].get_ylim();
	axarr[0].plot_date([EQtime, EQtime], [bottom, top], '--k');	
	axarr[0].set_ylabel("East (mm)");
	axarr[0].set_title("Detrended GPS Time Series")
	axarr[0].grid('on')

	for i in range(len(dataobj_list)):
		offset=spacing*i;
		ndata=dataobj_list[i].dN;
		ndata=[x + offset for x in ndata];

		l1 = axarr[1].plot_date(dataobj_list[i].dtarray,ndata,marker='+',markersize=2);
		line_color=l1[0].get_color()
		axarr[1].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,dataobj_list[i].name,fontsize=9,color=line_color);
	axarr[1].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20171020", "%Y%m%d"));
	axarr[1].set_ylim([-10,offset+10])
	bottom,top=axarr[1].get_ylim();
	axarr[1].plot_date([EQtime, EQtime], [bottom, top], '--k');	
	axarr[1].set_ylabel("North (mm)");
	axarr[1].set_title("Detrended GPS Time Series")
	axarr[1].grid('on')
	plt.savefig('Mend_Collective_TS_'+filename+'.jpg')	
	plt.close();


	return;




