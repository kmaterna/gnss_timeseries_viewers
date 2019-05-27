# The purpose of this script is to make a file that my GRACE code can use for one-time pre-computing of GRACE loading time series. 

import stations_within_radius
import gps_io_functions
import datetime as dt 

# Search for a box within northern California
# coord_box=[-125.5,-120,38,42];

coord_box=[-125.5,-119,36,46.5];
close_stations = stations_within_radius.get_stations_within_box(coord_box, network='pbo');
print(close_stations);

# Write them into a file. 
outfile=open("CA_OR.txt",'w');
for i in close_stations:
	[mydata] = gps_io_functions.read_pbo_pos_file("../../GPS_POS_DATA/PBO_Data/"+i+".pbo.final_nam08.pos");

	# Write into output file
	# Format: name, lon, lat, start, end

	starttime=dt.datetime.strftime(mydata.dtarray[0],"%Y");
	endtime=dt.datetime.strftime(mydata.dtarray[-1],"%Y");
	starttime=float(starttime);
	endtime=float(endtime)+1;
	print(endtime);
	if endtime>2019:
		endtime=2019;
	if starttime<2002:
		starttime=2002;
	outfile.write("%s %f %f %.2f %.2f\n" % (i, mydata.coords[0], mydata.coords[1], starttime, endtime) );


outfile.close();
