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
import single_station_tsplot
import stations_within_radius


def driver():
	[stations, distances, filenames, EQtime] = configure();
	dataobj_list = inputs(stations, filenames);

	[sorted_objects, sorted_distances] = compute(dataobj_list, distances);
	output_full_ts(sorted_objects, sorted_distances);

	return;

def configure():
	EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
	EQtime  = dt.datetime.strptime("20140310", "%Y%m%d");
	radius=150;  # km. 
	stations, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius);
	filenames=[];
	for station in stations:
		filenames.append("../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos");
	return [stations, distances, filenames, EQtime];

def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;


def compute(dataobj_list, distances):

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];

	new_data_obj=[];
	

	return [sorted_objects, sorted_distances];


def output_full_ts(dataobj_list, distances):

	plt.figure();
	[f,axarr]=plt.subplots(1,2,sharex=True,sharey=True)
	label_date="20171031";
	offset=0;
	for i in range(len(dataobj_list)):
		offset=16*i;
		edata=dataobj_list[i].dE;
		data_offset=[x + offset for x in edata];
		l1 = axarr[0].plot_date(dataobj_list[i].dtarray,edata,marker='+',markersize=2);
		line_color=l1[0].get_color()
		axarr[0].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,dataobj_list[i].name,fontsize=9,color=line_color);
	axarr[0].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20171020", "%Y%m%d"));
	axarr[0].set_ylim([-10,offset+10])
	axarr[0].set_ylabel("East (mm)");
	axarr[0].set_title("Detrended GPS Time Series")
	axarr[0].grid('on')

	for i in range(len(dataobj_list)):
		offset=16*i;
		ndata=dataobj_list[i].dN;
		data_offset=[x + offset for x in ndata];

		l1 = axarr[1].plot_date(dataobj_list[i].dtarray,ndata,marker='+',markersize=2);
		line_color=l1[0].get_color()
		axarr[1].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,dataobj_list[i].name,fontsize=9,color=line_color);
	axarr[1].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20171020", "%Y%m%d"));
	axarr[1].set_ylim([-10,offset+10])
	axarr[1].set_ylabel("North (mm)");
	axarr[1].set_title("Detrended GPS Time Series")
	axarr[1].grid('on')
	plt.savefig('Mend_Collective_TS.jpg')	
	plt.close();


	return;




