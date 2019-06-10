# simple script to print when stations have offsets larger than a certain value. 

import sys

inputfile=sys.argv[1];
cutoff=float(sys.argv[2]);

print("Reading for stations in %s above %fmm" % (inputfile, cutoff) );

ifile=open(inputfile,'r');
for line in ifile:
	temp=line.split();
	station=temp[0];
	east=float(temp[3]);
	north=float(temp[4]);

	if abs(east)>cutoff:
		print(station);
ifile.close();