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
import gps_ts_functions



def onset_time_map(name,lon,lat,data,size,earthquake_time,description,savename):
	
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
	map.drawcoastlines(linewidth=0.5);
	# map.fillcontinents(color='azure')
	map.drawmapboundary(fill_color='tan')
	# draw lat/lon grid lines every degree.
	map.drawmeridians(np.arange(-128,-120,1),labels=[1,0,0,1])
	map.drawparallels(np.arange(38,43,1),labels=[1,0,0,1])

	if description=='Up':
		multiplier=2;
	else:
		multiplier=10;

	# Plot dots
	for i in range(len(lon)):
		x,y=map(lon[i],lat[i])
		line_color=custom_cmap.to_rgba(data_points[i]);
		map.plot(x,y,marker='D',color=line_color,markersize=multiplier*size[i],linewidth=0.25);
		x,y=map(lon[i]+0.04, lat[i]-0.03);
		plt.text(x,y,name[i], fontsize=8);

	# Plot legend with earthquake time. 
	x,y=map(-124.8, 41.45)
	line_color=custom_cmap.to_rgba(0);
	map.plot(x,y,marker='D',color=line_color);
	x,y=map(-124.75, 41.41)
	plt.text(x,y,'Earthquake Time');

	custom_cmap.set_array(np.arange(vmin,vmax));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Inflection Time (Days since earthquake)');

	plt.title('GPS Onset Time near '+earthquake_time+', '+description)
	plt.savefig(savename);
	plt.close();
	return;



#  THE MAIN PROGRAM

earthquake_time="20140314"
# earthquake_time="20161208"
infile="Outputs/"+earthquake_time+"_inflections.txt"
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
	east_change.append(abs(float(temp[9])));
	north_change.append(abs(float(temp[10])));
	up_change.append(abs(float(temp[11])));

onset_time_map(name,lon,lat,east,east_change,earthquake_time,'East',earthquake_time+'_east.png');
onset_time_map(name,lon,lat,north,north_change,earthquake_time,'North',earthquake_time+'_north.png');
onset_time_map(name,lon,lat,up,up_change,earthquake_time,'Up',earthquake_time+'_up.png');


