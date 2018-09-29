import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from mpl_toolkits.basemap import Basemap


client_search="NCEDC"; # NCEDC has all the data; IRIS only has a little. 
fmin=0.01; fmax=0.2;
num_minutes=8;
starttime="2016-12-08T14:49:38.000"; # 2016 mendocino fault zone
# starttime="2016-12-08T17:38:46.000"; # solomon
# starttime="2014-03-10T05:17:54.000"; # 2014: no data
# starttime="2010-01-10T00:27:15.000"; # 2010. 
# starttime="2005-06-15T02:51:00.000"; # 2005
# starttime="2017-09-08T04:54:19.000"; # mexico
channel_type="HH*"
station_name="YBH";
network="BK";


print("working on %s" % station_name);
myclient=Client(client_search);
t = UTCDateTime(starttime);
print("Getting channel by searching %s" % (channel_type) );
print("myclient.get_waveforms(%s, %s, *, %s, %s, attach_response=True)\n" % (network, station_name, channel_type, starttime) );
try:
	st = myclient.get_waveforms(network, station_name, "*", channel_type, t, t + num_minutes * 60, attach_response=True);
except:
	print("could not find data. Exiting...");
	sys.exit(0);


st.remove_response(output="VEL");
st.filter("bandpass", freqmin=fmin, freqmax=fmax);
print(str(st));
print("\n");
print(str(st[0]));


transverse=st[0];

print(len(st[0].data))
print(len(st[1].data))
print(len(st[2].data))

plt.figure();
plt.plot(st[0].data,'b')
plt.plot(st[1].data,'r')
plt.plot(st[2].data,'k')
plt.savefig(station_name+'.eps');
plt.close();

print(np.nanmax(transverse.data));
print(np.nanmin(transverse.data));
print(transverse);

