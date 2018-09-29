# Plot all the available stations with their amplitude of shaking as their color. 

import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.imaging.beachball import beach
import sys
from mpl_toolkits.basemap import Basemap
import haversine



def run_map(eqcode):
	[lonbounds, latbounds, network_search, client_search, fmin, fmax, num_minutes, eqfile, ampfile]=configure(eqcode);
	[eqlon, eqlat, starttime, endtime]=inputs_one_earthquake(eqcode,eqfile);
	[namecodes, lons, lats, radials, transverses, verticals]=compute(lonbounds, latbounds, network_search, client_search, eqlon, eqlat, starttime, endtime, fmin, fmax, num_minutes);
	write_amplitudes(namecodes, lons, lats, radials, transverses, verticals, ampfile);
	amplitude_map(eqfile, ampfile, lonbounds, latbounds,eqcode,'verticals','Amplitude_maps/'+eqcode); # THE BASEMAP WITH STATIONS
	amplitude_map(eqfile, ampfile, lonbounds, latbounds,eqcode,'radials','Amplitude_maps/'+eqcode); # THE BASEMAP WITH STATIONS
	amplitude_map(eqfile, ampfile, lonbounds, latbounds,eqcode,'tangentials','Amplitude_maps/'+eqcode); # THE BASEMAP WITH STATIONS
	return;

def configure(eqcode):
	latbounds=[38.8,42];
	lonbounds=[-126.5,-121];
	eqfile="Amplitude_maps/eq_database.txt";
	ampfile="Amplitude_maps/outputs_"+eqcode+".txt";
	client_search="NCEDC"; # NCEDC has all the data; IRIS only has a little. 
	network_search="BK,NC,PB"; # LOTS OF DATA
	# network_search="NC";
	fmin=0.005; fmax=0.05; # surface waves
	# fmin=0.1; fmax=2.0;  # p waves 
	num_minutes=6;	  # 7 for local events. 
	if eqcode=="2017" or eqcode=="2018":
		num_minutes=60; # 20 for teleseismic event. 
	return [lonbounds, latbounds, network_search, client_search, fmin, fmax, num_minutes, eqfile,ampfile];




def inputs_one_earthquake(eqcode,eqfile):
	# Reading in reference info about the earthquake. 
	ifile=open(eqfile,'r');
	for line in ifile:
		temp=line.split();
		if eqcode in temp[0]:
			eqlon=float(temp[2]);
			eqlat=float(temp[3]);
			starttime=temp[4];
			endtime=temp[5];
	ifile.close();
	return [eqlon, eqlat, starttime, endtime];

def inputs_all_earthquakes(eqfile):
	datelabelarray=[]; magarray=[]; eqlonarray=[]; eqlatarray=[]; mtarray=[];
	ifile=open(eqfile);
	for line in ifile:
		temp=line.split();
		datelabelarray.append(temp[0]);
		magarray.append(temp[1]);
		eqlonarray.append(float(temp[2]));
		eqlatarray.append(float(temp[3]));
		mtarray.append([float(temp[6]), float(temp[7]), float(temp[8])]);
	return [datelabelarray, magarray, eqlonarray, eqlatarray, mtarray];

def inputs_amplitudes(ampfile):
	name=[]; lon=[]; lat=[]; radials=[]; transverses=[]; verticals=[];
	ifile=open(ampfile,'r');
	for line in ifile:
		temp=line.split();
		name.append(temp[0]);
		lon.append(float(temp[1]));
		lat.append(float(temp[2]));
		radials.append(float(temp[3]));
		transverses.append(float(temp[4]));
		verticals.append(float(temp[5]));
	ifile.close();
	return [name, lon, lat, radials, transverses, verticals];

def write_amplitudes(namecodes, lons, lats, radials, transverses, verticals, ampfile):
	ofile=open(ampfile,'w');
	for i in range(len(lons)):
		ofile.write('%s %f %f %f %f %f\n' % (namecodes[i], lons[i],lats[i],radials[i], transverses[i], verticals[i] ) );
	ofile.close();
	return

def make_plot(st, channel_id, station):
	# Plot the waveforms. 
	plt.figure();
	if len(st)==1:
		plt.plot(st[0],'k');
	elif len(st)==2:
		plt.plot(st[0],'k');
	elif len(st)==3:
		ymax=np.max([np.max(st[0].data),np.max(st[1].data)]);
		f, axarr=plt.subplots(3,1,sharex=True);
		axarr[0].plot(st[0],'b');
		axarr[0].text(0,0,'transverse',color='blue');
		axarr[0].set_ylim([-ymax,ymax]);
		axarr[1].plot(st[1],'r');
		axarr[1].text(0,0,'radial',color='red');
		axarr[1].set_ylim([-ymax,ymax]);
		axarr[2].plot(st[2],'k');
		axarr[2].text(0,0,'vertical',color='black');
		axarr[2].set_ylim([-ymax,ymax]);
	plt.savefig('Amplitude_maps/'+channel_id+'_'+station+'.png');
	plt.close();
	return;



def compute(lonbounds, latbounds, network_search, client_search, eqlon, eqlat, starttime, endtime, fmin, fmax, num_minutes):
	myclient=Client(client_search);
	waveform_client=Client("NCEDC");  # this client has more waveforms than the IRIS client. 
	network_split=network_search.split(',');  # [BK,NC,etc]. 

	inventory = myclient.get_stations(network=network_search, starttime=starttime, endtime=endtime, 
	minlatitude=latbounds[0], maxlatitude=latbounds[1], minlongitude=lonbounds[0], maxlongitude=lonbounds[1], 
	level='response', channel='?H?');  # ?H? means broadband. 
	# len(inventory) is the number of networks. inventory[0] is associated with the first network code.

	eqtuple=(eqlat,eqlon); # this is the required order for azimuth math. 
	

	namecodes=[]; lons=[]; lats=[]; radials=[]; transverses=[]; verticals=[];

	for i in range(len(inventory)):
		print(inventory[i]);

	for i in range(len(inventory)):  # for each network
		for j in range(len(inventory[i])):  # flipping through stations
			new_st=[];
			newcount=0;


			mystation=inventory[i][j];  # a station object. 
			station_dictionary=mystation.get_contents();
			station_name=station_dictionary['stations'][0].split(' ')[0];
			# if station_name != 'GASB' and station_name != 'KHMB' and station_name != 'SUTB' and station_name != 'HATC':
				# continue;		# GET RID OF THIS LATER. 	

			print("working on %s" % station_name);

			all_channels=station_dictionary['channels'];
			print(all_channels);
			channel0=all_channels[0]; # the first channel in the list (not very specific).
			channel_id=network_split[i]+'.'+channel0;
			channel_type=channel_id.split('.')[-1][0:2]+"*";  # something like "BH*" or "HH*" depending on what the station has. 

			# Parameters
			coords=inventory.get_coordinates(channel_id);
			pt_station=(coords['latitude'],coords['longitude']);
			back_az=haversine.calculate_initial_compass_bearing(pt_station, eqtuple);
			print(pt_station);
			print(eqtuple);
			print(back_az);

			t = UTCDateTime(starttime);	# Parameters

			# Get some waveforms
			print("getting waveforms:");
			print("waveform_client.get_waveforms(%s, %s, *, %s, %s, attach_response=True)\n" % (network_split[i], station_name, channel_type, starttime) );
			try:
				st = waveform_client.get_waveforms(network_split[i], station_name, "*", channel_type, t, t + num_minutes * 60, attach_response=True);
			except:
				print("Could not find data for station %s. Skipping. " % (station_name) );
				continue;

			if station_name=='JCC':
				print("waveform_client.get_waveforms(%s, %s, *, %s, %s, attach_response=True)\n" % (network_split[i], station_name, "HH*", starttime) );
				try:
					st = waveform_client.get_waveforms(network_split[i], station_name, "*", "HH*", t, t + num_minutes * 60, attach_response=True);
				except:
					print("Could not find data for station %s. Skipping. " % (station_name) );
					continue;

			obtained_channel=[];
			for k in range(len(st)):
				obtained_channel.append(str(st[k]).split()[0]);  # appends things like 'BK.YBH.00.BHE'. 
			print("Got channels:");
			print(obtained_channel);

			# Get the data from the station, and record its amplitude. 
			# This is the meat of the program. 

			# Processing once we have the stream. 
			# Vertical only stations
			if len(st)==1 and obtained_channel[0][-1]=='Z':
				print("Only vertical traces");
				st.remove_response(output="VEL");
				st.filter("bandpass", freqmin=fmin, freqmax=fmax);
				vertical=st[0];
				vertical_number=np.nanmax(vertical);
				radial_number=np.nan;
				transverse_number=np.nan;
				# make_plot(st, channel_id, station_name);

			elif len(st)==2:
				if obtained_channel[0][-1]=='Z':
					vertical=st[0];
				elif obtained_channel[1][-1]=='Z':
					vertical=st[1];
				else:
					print("cannot find vertical; continuing");
					continue;
				print("Only using one vertical trace");
				st.remove_response(output="VEL");
				st.filter("bandpass", freqmin=fmin, freqmax=fmax);
				vertical=st[0];
				vertical_number=np.nanmax(vertical);
				radial_number=np.nan;
				transverse_number=np.nan;
				# make_plot(st, channel_id, station_name);
			
			# Three component stations.
			elif len(st)==3:
				# We had a weird situation where JCC.BH* didn't have the same length of data. Trying a hard-coded fix with HH*. 
				if len(st[0].data)!=len(st[1].data) or len(st[1].data)!=len(st[2].data) or len(st[0].data)!=len(st[2].data): 
					print("Bad data.");
					try:
						st = waveform_client.get_waveforms(network_split[i], station_name, "*", "HH*", t, t + num_minutes * 60, attach_response=True);
					except:
						print("Could not find HH* data for station %s. Skipping. " % (station_name) );
						continue;

				# Remove the response and rotate to transverse. 
				st.remove_response(output="VEL");
				st.rotate('NE->RT',back_azimuth=back_az);
				st.filter("bandpass", freqmin=fmin, freqmax=fmax);
				make_plot(st,channel_id,station_name);
				
				transverse=st[0];
				radial=st[1];
				vertical=st[2];
				radial_number=np.nanmax(radial);
				transverse_number=np.nanmax(transverse);
				vertical_number=np.nanmax(vertical);

			elif len(st)>3:
				print("Have %d traces. Removing some." % len(st));
				for k in range(len(st)):
					if '.00.' in obtained_channel[k]:
						newcount=newcount+1;
				print("Found %d trances with .00." % newcount);	

				if station_name=='JCC': # JCC is so messy!  
					if len(st)==6:
						new_st.append(st[0]);
						new_st.append(st[2]);
						new_st.append(st[4]);
					# 2014 has 4 traces with different time spans. This isn't so great. 
					else:
						print("Don't know what to do with JCC.");
						continue;
				elif newcount==3:
					for k in range(len(st)):
						if '.00.' in obtained_channel[k]:
							new_st.append(st[k]);
				else:
					print("Didn't find 3 .00. channels. Don't know what to do.");
					continue;

				# Remove the response and rotate to transverse. 
				st.remove_response(output="VEL");
				st.rotate('NE->RT',back_azimuth=back_az);
				st.filter("bandpass", freqmin=fmin, freqmax=fmax);
				# make_plot(st,channel_id,station_name);
				
				transverse=st[0];
				radial=st[1];
				vertical=st[2];
				radial_number=np.nanmax(radial);
				transverse_number=np.nanmax(transverse);
				vertical_number=np.nanmax(vertical);

			else:
				print("stream is %d traces long... not sure what to do" % (len(st) ) );
				print(st);
				continue;

			# Saving data for the plotting function. 
			lons.append(coords['longitude']);
			lats.append(coords['latitude']);
			namecodes.append(station_name);
			transverses.append(transverse_number);
			radials.append(radial_number);
			verticals.append(vertical_number);

	return [namecodes, lons, lats, radials, transverses, verticals];






def amplitude_map(eqfile, ampfile, lonbounds,latbounds,eqcode, component, savename):
	
	[datelabelarray, magarray, eqlonarray, eqlatarray, mtarray]=inputs_all_earthquakes(eqfile);
	[name, lon, lat, radials, transverses, verticals]=inputs_amplitudes(ampfile);

	if component=='verticals':
		plotting_array = verticals;
	elif component=='radials':
		plotting_array = radials;
	elif component=='tangentials':
		plotting_array = transverses;

	# Restricted colorscale
	vmin=np.nanmin(plotting_array);
	vmax=np.nanmax(plotting_array);
	vmin=0;
	vmax=2;

	color_boundary_object=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet');

	# draw coastlines, country boundaries, fill continents.
	plt.figure(figsize=(10,10));
	map=Basemap(projection='merc',llcrnrlat=latbounds[0],llcrnrlon=lonbounds[0], urcrnrlat=latbounds[1], urcrnrlon=lonbounds[1],resolution='i');
	map.drawcoastlines(linewidth=0.5);
	map.fillcontinents(color='tan',lake_color='azure')
	map.drawmapboundary(fill_color='azure')
	# draw lat/lon grid lines every degree.
	map.drawmeridians(np.arange(lonbounds[0],lonbounds[1],1),labels=[1,0,0,1]);
	map.drawparallels(np.arange(latbounds[0],latbounds[1],1),labels=[1,0,0,1]);

	# Plot dots colored by amplitude
	for i in range(len(lon)):
		if ~np.isnan(plotting_array[i]):
			x,y=map(lon[i],lat[i]);
			line_color=custom_cmap.to_rgba(plotting_array[i]);
			map.plot(x,y,marker='D',color=line_color,markeredgecolor='black',markeredgewidth=0.3); # markersize=multiplier*size[i]
			x,y=map(lon[i]+0.04, lat[i]-0.03);
			plt.text(x,y,name[i], fontsize=10);

	# Plot earthquakes
	for i in range(len(datelabelarray)):
		x,y=map(eqlonarray[i],eqlatarray[i]);
		b = beach(mtarray[i], xy=(x, y), width=30000, linewidth=1)
		plt.gca().add_collection(b);
		x,y=map(eqlonarray[i]-0.1, eqlatarray[i]-0.21);
		plt.text(x,y,magarray[i], fontsize=10);

	custom_cmap.set_array(np.arange(vmin,vmax));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Amplitude',fontsize=14);
	cb.ax.tick_params(labelsize=14)

	plt.title('Surface wave amplitude for '+component+' during '+eqcode,fontsize=16)
	plt.savefig(savename+'_'+component+'.png');
	plt.close();
	
	return;






if __name__=="__main__":
	# run_map("2005");
	# run_map("2010");
	run_map("2014");
	# run_map("2016");
	# run_map("2017");
	# run_map("2018");

