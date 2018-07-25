# Plotting the results

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime as dt 
import matplotlib
import matplotlib.cm as cm
import gps_ts_functions



def onset_time_map(lon,lat,data,earthquake_time,description,savename):
	
	# Earthquake time
	eq_dt = dt.datetime.strptime(earthquake_time,"%Y%m%d");

	# Make a color scale based on the time. 
	data_points = [(i-eq_dt).days for i in data];
	
	# Restricted colorscale
	threshold=180;  # days
	vmin=-threshold;
	vmax=threshold;

	color_boundary_object=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# draw coastlines, country boundaries, fill continents.
	map=Basemap(projection='merc',llcrnrlat=39,llcrnrlon=-125, urcrnrlat=41.5, urcrnrlon=-122,resolution='i');
	map.drawcoastlines(linewidth=0.25);
	map.fillcontinents(color='coral',lake_color='white')
	map.drawmapboundary(fill_color='white')
	# draw lat/lon grid lines every degree.
	map.drawmeridians(np.arange(-128,-120,1),labels=[1,0,0,1])
	map.drawparallels(np.arange(38,43,1),labels=[1,0,0,1])

	# Plot dots
	for i in range(len(lon)):
		x,y=map(lon[i],lat[i])
		line_color=custom_cmap.to_rgba(data_points[i]);
		map.plot(x,y,marker='D',color=line_color);

	# Plot legend with earthquake time. 
	x,y=map(-124.8, 41.45)
	line_color=custom_cmap.to_rgba(0);
	map.plot(x,y,marker='D',color=line_color);
	x,y=map(-124.75, 41.41)
	plt.text(x,y,'Earthquake Time');

	custom_cmap.set_array(np.arange(vmin,vmax));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Inflection Time (Days since earthquake)');

	plt.title('GPS Onset Time, '+description)
	plt.savefig(savename);
	plt.close();
	return;



#  THE MAIN PROGRAM

earthquake_time="20140314"
# earthquake_time="20161208"
infile="Outputs/"+earthquake_time+"_inflections.txt"
lat=[]; lon=[]; east=[]; north=[]; up=[];

ifile=open(infile,'r');
for line in ifile:
	temp=line.split();
	lon.append(float(temp[0]));
	lat.append(float(temp[1]));
	east.append(dt.datetime.strptime(temp[2],"%Y-%m-%d"));
	north.append(dt.datetime.strptime(temp[4],"%Y-%m-%d"));
	up.append(dt.datetime.strptime(temp[6],"%Y-%m-%d"));


onset_time_map(lon,lat,east,earthquake_time,'East',earthquake_time+'_east.png');
onset_time_map(lon,lat,north,earthquake_time,'North',earthquake_time+'_north.png');
onset_time_map(lon,lat,up,earthquake_time,'Up',earthquake_time+'_up.png');


