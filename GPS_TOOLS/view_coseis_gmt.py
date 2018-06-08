# Viewing offsets from a given earthquake
import numpy as np 
import glob as glob
import subprocess


def driver():
	filename=config();
	[lon, lat, dE, dN, sE, sN, dU, sU, site, eqlon, eqlat]=inputs(filename);
	outputs(lon, lat, dE, dN, sE, sN, dU, sU, site, eqlon, eqlat);
	return;

def config():
	events_dir="../../GPS_POS_DATA/Event_Files/"
	event = '140310'  # format: March 10, 2014  #THIS IS WHERE YOU DEFINE THE EARTHQUAKE YOU WANT
	# event = '100110'
	event = '050615'
	filenames=glob.glob(events_dir+'pbo_'+event+'*_coseis_kalts.evt');
	if len(filenames)!=1: 
		print("Error: The wrong number of files has been detected.")
		print(filenames);
	return filenames[0];

def inputs(filename):
	lon=[]; lat=[]; dE=[]; dN=[]; sE=[]; sN=[]; dU=[]; sU=[]; site=[]; 
	infile=open(filename,'r');
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
	return [lon, lat, dE, dN, sE, sN, dU, sU, site, eqlon, eqlat];


def outputs(lon, lat, dE, dN, sE, sN, dU, sU, site, eqlon, eqlat):
	#### Will eventually make these parameters. 
	lonW=-126.7
	lonE=-122.0
	latS=39.5
	latN=42.0
	####
	ofile=open('coseis_vectors.txt','w');
	for i in range(len(lon)):
		ofile.write(str(lon[i])+" "+str(lat[i])+" "+str(dE[i])+" "+str(dN[i])+" "+str(sE[i])+" "+str(sN[i])+" 0 "+site[i]+" \n");
	ofile.close();

	subprocess.call(['./coseis_map_gps.gmt',str(lonW),str(lonE),str(latS),str(latN),str(eqlon),str(eqlat)],shell=False);
	return;

if __name__=="__main__":
	driver();

