# Viewing coseismic offsets from a given earthquake. 
# Updated to view 3 earthquakes at once. 
# MUST VIEW ALL 3 EARTHQUAKES. 

import numpy as np 
import glob as glob
import subprocess


def driver():
	[filenames, events] =config();
	[lons, lats, dEs, dNs, sEs, sNs, dUs, sUs, sites, eqlons, eqlats]=inputs(filenames);
	outputs(lons, lats, dEs, dNs, sEs, sNs, dUs, sUs, sites, eqlons, eqlats, events);
	return;

def config():
	events_dir="../../GPS_POS_DATA/Event_Files/"
	events=['050615','100110','140310']; # format: YYMMDD
	filenames=[];
	for item in events:
		ifilename=glob.glob(events_dir+'pbo_'+item+'*_coseis_kalts.evt');
		filenames.append(ifilename[0]);
	if len(filenames)==0: 
		print("Error: No files have been detected.");
	return [filenames, events];

def inputs(filenames):  # can take an array of filenames and read them into a list of lists. 
	print(filenames);
	lons =[]; lats=[]; dEs=[]; dNs=[]; sEs=[]; sNs=[]; dUs=[]; sUs=[]; sites=[]; eqlons=[]; eqlats=[]; 
	for i in range(len(filenames)):
		lon=[]; lat=[]; dE=[]; dN=[]; sE=[]; sN=[]; dU=[]; sU=[]; site=[]; 
		infile=open(filenames[i],'r');
		read_data=0;
		for line in infile:
			temp=line.split();
			if "mm" in temp and "deg" in temp:
				read_data=1;
				continue;
			else:
				if "(lat/long)" in temp:
					eqlat=float(temp[4]);
					eqlon=float(temp[5])-360;
			if read_data==1:
				lon.append(float(temp[0])-360);
				lat.append(float(temp[1]));
				dE.append(float(temp[2]));
				dN.append(float(temp[3]));
				sE.append(float(temp[4]));
				sN.append(float(temp[5]));
				dU.append(float(temp[7]));
				sU.append(float(temp[8]));
				full_name=temp[9];
				site.append(full_name[0:4]);
		lons.append(lon); lats.append(lat); dEs.append(dE); dNs.append(dN); sEs.append(sE); sNs.append(sN); dUs.append(dU); sUs.append(sU); 
		sites.append(site); eqlons.append(eqlon); eqlats.append(eqlat);
	return [lons, lats, dEs, dNs, sEs, sNs, dUs, sUs, sites, eqlons, eqlats];


def outputs(lon, lat, dE, dN, sE, sN, dU, sU, site, eqlon, eqlat, events):
	#### Will eventually make these parameters. 
	lonW=-126.7
	lonE=-122.5
	latS=39.5
	latN=42.0
	####
	event_files=['','',''];
	for i in range(len(events)):
		ofile=open('coseis_vectors_'+events[i]+'.txt','w');
		event_files[i]='coseis_vectors_'+events[i]+'.txt';
		for j in range(len(lon[i])):
			ofile.write(str(lon[i][j])+" "+str(lat[i][j])+" "+str(dE[i][j])+" "+str(dN[i][j])+" "+str(sE[i][j])+" "+str(sN[i][j])+" 0 "+site[i][j]+" \n");
		ofile.close();
	subprocess.call(['./coseis_map_gps.gmt',str(lonW),str(lonE),str(latS),str(latN),event_files[0],event_files[1], event_files[2]],shell=False);
	return;

if __name__=="__main__":
	driver();

