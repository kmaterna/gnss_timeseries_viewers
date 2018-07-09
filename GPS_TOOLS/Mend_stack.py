# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and eventually trend from the data. 
# Step 4: Plot in order of increasing latitude, colored by how close they are to the earthquake. 

# Reference: 
# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import collections
import datetime as dt 
from scipy import signal
import gps_io_functions
import gps_ts_functions
import stations_within_radius


def driver():
	[stations, distances, filenames, EQtime, earthquakes_dir, offsets_dir] = configure();
	dataobj_list = inputs(stations, filenames);
	[sorted_objects, noeq_objects, sorted_distances, east_slope_obj] = compute(dataobj_list, distances, earthquakes_dir, offsets_dir, EQtime);
	
	output_full_ts(sorted_objects, sorted_distances, EQtime, "detrended", east_slope_obj);
	output_full_ts(noeq_objects, sorted_distances, EQtime, "noeq_noseasons", east_slope_obj);
	return;


def configure():
	EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
	EQtime  = dt.datetime.strptime("20140310", "%Y%m%d");
	earthquakes_dir = "../GPS_POS_DATA/Event_Files/";
	offsets_dir = "../GPS_POS_DATA/Offsets/";
	radius=120;  # km. 
	stations, distances = stations_within_radius.get_stations_within_radius(EQcoords, radius);
	filenames=[];
	for station in stations:
		filenames.append("../GPS_POS_DATA/PBO_stations/"+station+".pbo.final_nam08.pos");
	return [stations, distances, filenames, EQtime, earthquakes_dir, offsets_dir];

def inputs(stations, filenames):
	dataobj_list=[];
	for item in filenames:
		[myData]=gps_io_functions.read_pbo_pos_file(item);
		dataobj_list.append(myData);
	return dataobj_list;


def compute(dataobj_list, distances, earthquakes_dir, offsets_dir, EQtime):

	latitudes_list=[i.coords[1] for i in dataobj_list];
	sorted_objects = [x for _,x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data. 
	sorted_distances = [x for _,x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

	# Detrended objects
	detrended_objects=[];
	for i in range(len(dataobj_list)):
		newobj=gps_ts_functions.detrend_data(sorted_objects[i]);
		detrended_objects.append(newobj);

	# Objects with no earthquakes or seasonals
	noeq_objects = [];
	east_slope_obj=[];
	for i in range(len(dataobj_list)):
		# Remove the earthquakes
		newobj=gps_ts_functions.remove_offsets(sorted_objects[i],offsets_dir);
		newobj=gps_ts_functions.remove_earthquakes(newobj,earthquakes_dir);
		newobj=gps_ts_functions.remove_annual_semiannual(newobj);

		# Get the pre-event and post-event velocities (earthquakes removed)
		[east_slope_before, north_slope_before, vert_slope_before]=gps_ts_functions.get_slope(newobj,endtime=EQtime);
		[east_slope_after, north_slope_after, vert_slope_after]=gps_ts_functions.get_slope(newobj,starttime=EQtime);
		east_slope_after=np.round(east_slope_after,decimals=1);
		east_slope_before=np.round(east_slope_before,decimals=1);
		east_slope_obj.append([east_slope_before, east_slope_after]);

		# Remove trend and then save this new object
		newobj=gps_ts_functions.detrend_data(newobj);
		noeq_objects.append(newobj);


	return [detrended_objects, noeq_objects, sorted_distances, east_slope_obj];




def output_full_ts(dataobj_list, distances, EQtime, filename, east_slope_obj):

	plt.figure();
	[f,axarr]=plt.subplots(1,2,sharex=True,sharey=True)
	label_date="20181031";
	EQ1time = dt.datetime.strptime("20050615", "%Y%m%d");
	EQ2time = dt.datetime.strptime("20100110", "%Y%m%d");
	offset=0;
	spacing=10;
	closest_station=70;  # km from event
	farthest_station=120; # km from event
	color_boundary_object=matplotlib.colors.Normalize(vmin=closest_station,vmax=farthest_station, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	for i in range(len(dataobj_list)):
		offset=spacing*i;
		edata=dataobj_list[i].dE;
		edata=[x + offset for x in edata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = axarr[0].plot_date(dataobj_list[i].dtarray,edata,marker='+',markersize=2,color=line_color);
		axarr[0].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,dataobj_list[i].name,fontsize=9,color=line_color);
		# axarr[0].text(dt.datetime.strptime("20050301", "%Y%m%d"),offset,east_slope_obj[i][0],fontsize=9,color='k');
		# axarr[0].text(EQtime,offset,east_slope_obj[i][1],fontsize=9,color='k');
	axarr[0].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20181020", "%Y%m%d"));
	axarr[0].set_ylim([-10,offset+10])
	bottom,top=axarr[0].get_ylim();
	axarr[0].plot_date([EQtime, EQtime], [bottom, top], '--k');	
	axarr[0].plot_date([EQ1time, EQ1time], [bottom, top], '--k');	
	axarr[0].plot_date([EQ2time, EQ2time], [bottom, top], '--k');	
	axarr[0].set_ylabel("East (mm)");
	axarr[0].set_title("Detrended GPS Time Series")
	axarr[0].grid('on')

	for i in range(len(dataobj_list)):
		offset=spacing*i;
		ndata=dataobj_list[i].dN;
		ndata=[x + offset for x in ndata];
		line_color=custom_cmap.to_rgba(distances[i]);
		l1 = axarr[1].plot_date(dataobj_list[i].dtarray,ndata,marker='+',markersize=2, color=line_color);
		#axarr[1].text(dt.datetime.strptime(label_date, "%Y%m%d"),offset,dataobj_list[i].name,fontsize=9,color=line_color);
	axarr[1].set_xlim(dt.datetime.strptime("20050101", "%Y%m%d"),dt.datetime.strptime("20181020", "%Y%m%d"));
	axarr[1].set_ylim([-10,offset+10])
	bottom,top=axarr[1].get_ylim();
	axarr[1].plot_date([EQtime, EQtime], [bottom, top], '--k');	
	axarr[1].plot_date([EQ1time, EQ1time], [bottom, top], '--k');	
	axarr[1].plot_date([EQ2time, EQ2time], [bottom, top], '--k');		
	axarr[1].set_ylabel("North (mm)");
	axarr[1].set_title("Detrended GPS Time Series")
	axarr[1].grid('on')
	custom_cmap.set_array(range(closest_station,farthest_station))
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilometers from 2014 Earthquake');
	plt.savefig('Mend_Collective_TS_'+filename+'.jpg')	
	plt.close();


	return;




