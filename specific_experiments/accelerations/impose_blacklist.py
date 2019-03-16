# Python script to remove stations from pre-computed files if they have just been added to the blacklist
# Note that when you do this, it's deleting data
# It's going to be annoying to get that data back. 

import numpy as np
import glob
import sys
from subprocess import call


def remove_blacklist(directory, blacklist_file):
	bad_stations = read_blacklist(blacklist_file);
	file_list = glob.glob(directory+"*.txt");
	for datafile in file_list:
		print("Removing stations from %s " % datafile);
		print(bad_stations);
		remove_stations(datafile, bad_stations);
	return;

def read_blacklist(ifile):
	infile = open(ifile,'r');
	bad_stations=[];
	for line in infile:
		temp=line.split();
		bad_stations.append(temp[0]);
	infile.close();
	return bad_stations;

def remove_stations(datafile, bad_stations):
	with open(datafile,'r') as content_file:
		content=content_file.read();
	
	ofile=open(datafile,'w');
	for line in content.split('\n'):
		if len(line.split())==0:
			continue;
		if line.split()[0]=="#":
			ofile.write(line+"\n");
		else:
			if line.split()[-1] not in bad_stations:
				ofile.write(line+"\n");
			else:
				print("Removing "+line.split()[-1]);
	ofile.close();
	return;



if __name__=="__main__":
	directory = "unr_lssq_ITRF/"
	blacklist = "../../GPS_POS_DATA/Blacklist/blacklist.txt"
	remove_blacklist(directory, blacklist);
