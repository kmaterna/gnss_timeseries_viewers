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
	plt.figure(figsize=(8,8));
	map=Basemap(projection='merc',llcrnrlat=39,llcrnrlon=-125, urcrnrlat=41.5, urcrnrlon=-122,resolution='i');
	map.drawcoastlines(linewidth=0.5);
	map.fillcontinents(color='tan',lake_color='azure')
	map.drawmapboundary(fill_color='azure')
	# draw lat/lon grid lines every degree.
	map.drawmeridians(np.arange(-128,-120,1),labels=[1,0,0,1])
	map.drawparallels(np.arange(38,43,1),labels=[1,0,0,1])

	if description=='Vertical':
		multiplier=2;
	else:
		multiplier=10;

	# Plot dots
	for i in range(len(lon)):
		x,y=map(lon[i],lat[i])
		line_color=custom_cmap.to_rgba(data_points[i]);
		map.plot(x,y,marker='D',color=line_color,markersize=multiplier*size[i],markeredgecolor='black',markeredgewidth=0.3);
		x,y=map(lon[i]+0.04, lat[i]-0.03);
		plt.text(x,y,name[i], fontsize=10);

	# Plot legend with earthquake time. 
	x,y=map(-124.8, 41.45)
	line_color=custom_cmap.to_rgba(0);
	map.plot(x,y,marker='D',color=line_color,markeredgecolor='black',markeredgewidth=0.3);
	x,y=map(-124.75, 41.41)
	plt.text(x,y,'Earthquake Time',fontsize=14);

	custom_cmap.set_array(np.arange(vmin,vmax));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Onset Time (Days since earthquake)',fontsize=14);
	cb.ax.tick_params(labelsize=14)

	plt.title(description+' GPS Onset Time near '+earthquake_time,fontsize=16)
	plt.savefig(savename);
	plt.close();
	return;


def onset_time_map_GMT(name,lon,lat,east,east_change,north, north_change, vert, vert_change, earthquake_time):
	# Earthquake time
	eq_dt = dt.datetime.strptime(earthquake_time,"%Y%m%d");
	max_size=3.5;
	outfile_name=earthquake_time+".txt";
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


	call(['./timing_map_gps.gmt',outfile_name,'-125', '-118', '36.5', '42.0',earthquake_time+"_norcal",earthquake_time],shell=False);  # northern California
	call(['./timing_map_gps.gmt',outfile_name,'-125', '-110', '32.5', '48.5',earthquake_time+"_WUS",earthquake_time],shell=False);   # WUS
	call(['./timing_map_gps.gmt',outfile_name,'-122', '-115', '32.5', '37.5',earthquake_time+"_socal", earthquake_time],shell=False);   # SoCal
	return;



if __name__=="__main__":

	#  THE MAIN PROGRAM
	# earthquake_time="20140310"  #
	earthquake_time="20161208";
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
		east_change.append(float(temp[9]));
		north_change.append(float(temp[10]));
		up_change.append(float(temp[11]));

	# Make GMT maps. 
	onset_time_map_GMT(name, lon, lat, east, east_change, north, north_change, up, up_change, earthquake_time);

