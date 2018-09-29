# A little script to download data from a station, and plot it for all 4 events. 

#obspy client for get_events(), get_stations(), get_waveforms(). 
#obspy waveform fetch
#obspy mass downloader * can be used for getting data by events. 
#obspy remove_response() 


import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import sys
import haversine


# Reference information. 
client = Client("NCEDC")
station="KHMB"
# station="KMPB"
pt_station=(40.2022, -123.7086); # kmr
num_minutes=7.0;
fmin=0.1;
fmax=1.0;
plotname=station+'_'+str(fmin)+'_'+str(fmax)+'.png'


# The event information. 
event_tuples=[]; time_list=[]; name_list=[];
# 2005 event
event_tuples.append((41.292, -125.952)); 
time_list.append("2005-06-15T02:51:00.000"); 
name_list.append("2005 M7.2");

# 2010 event
event_tuples.append((40.652, -125.692));
time_list.append("2010-01-10T00:27:15.000");
name_list.append("2010 M6.5");

# 2014 event
event_tuples.append((40.829, -125.134)); 
time_list.append("2014-03-10T05:17:54.000"); 
name_list.append("2014 M6.8");

# 2016 event
event_tuples.append((40.453, -126.194)); 
time_list.append("2016-12-08T14:49:38.000"); 
name_list.append("2016 M6.6");





fig,axs=plt.subplots(4,1,sharex=True,figsize=(10,10));
fig.subplots_adjust(hspace=0);

# The waveform plots
xarray=[]; radial=[]; transverse=[]; vertical=[];
for i in range(len(event_tuples)):
	print("working on %s" % name_list[i])
	back_az=haversine.calculate_initial_compass_bearing(pt_station, event_tuples[i]);
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


