# How does a regional M6.6 compare to Tohoku? 

# Main tasks: 
# Read catalog
# Collect data
# Write data
# 


import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import sys
import haversine

def configure():
	input_file="major_eq_cat.txt";
	client = Client("NCEDC")
	# station="KHMB"
	station="KMPB"
	# station="KMR"
	station_coords=(40.2022, -123.7086); # kmr
	num_minutes=90.0;
	fmin=0.5;
	fmax=3.5;	
	return [client, station, station_coords, num_minutes, fmin, fmax, input_file]; 

def inputs(input_file):
	ifile=open(input_file,'r');
	event_tuples=[]; time_list=[]; name_list=[];
	for line in ifile:
		temp=line.split();
		lon=float(temp[0]);
		lat=float(temp[1]);
		event_tuples.append((lat, lon));
		time_list.append(temp[2]);
		name_list.append(temp[4:]);
		if len(name_list)==1:
			break;
	ifile.close();
	return [event_tuples, time_list, name_list];

def plot(client, station, station_coords, num_minutes, fmin, fmax, event_tuples, time_list, name_list):
	plotname=station+'_'+str(fmin)+'_'+str(fmax)+'_data.png'
	fig,axs=plt.subplots(4,1,sharex=True,figsize=(10,10));
	fig.subplots_adjust(hspace=0);

	# The waveform plots
	xarray=[]; radial=[]; transverse=[]; vertical=[];
	for i in range(len(event_tuples)):
		print("working on %s" % name_list[i])
		back_az=haversine.calculate_initial_compass_bearing(station_coords, event_tuples[i]);
		t = UTCDateTime(time_list[i])
		st = client.get_waveforms("NC", station, "*", "HH*", t, t + num_minutes * 60, attach_response=True);
		st.remove_response(output="VEL");
		st.rotate('NE->RT',back_azimuth=back_az);
		tr=st[0];  # the first trace
		df=tr.stats.sampling_rate;
		npts=tr.stats.npts;
		st.filter("bandpass", freqmin=fmin, freqmax=fmax);
		transverse.append(st[0]);
		radial.append(st[1]);
		vertical.append(st[2]);


	tempmax=[]; 
	for i in range(len(event_tuples)):
		tempmax.append(max(radial[i])+max(transverse[i])+max(vertical[i]))
	realmax=max(tempmax);
	x=np.arange(0,npts/df,1.0/df);

	for i in range(len(event_tuples)):
		newradial = [x+realmax for x in radial[i]];
		newvertical = [x+realmax*1.5 for x in vertical[i]];
		h1=axs[i].plot(x,newvertical,label='vertical');
		h2=axs[i].plot(x,newradial,label='radial');
		h3=axs[i].plot(x,transverse[i],label='transverse')
		axs[i].set_ylim(-realmax, 2.2*realmax);
		axs[i].set_ylabel('Vel (m/s)');
		axs[i].text(0,1.7*realmax,name_list[i],fontsize=16);
		if i==0:
			axs[i].legend(loc=1,fontsize=14);
	axs[-1].set_xlabel('Time (s)');
	axs[0].set_title(station+' Records at %.2f-%.2f Hz' % (fmin, fmax) ,fontsize=20);
	plt.savefig(plotname);

	return;

if __name__=="__main__":
	[client, station, station_coords, num_minutes, fmin, fmax, input_file] = configure();
	[event_tuples, time_list, name_list] = inputs(input_file);
	plot(client, station, station_coords, num_minutes, fmin, fmax, event_tuples, time_list, name_list);

