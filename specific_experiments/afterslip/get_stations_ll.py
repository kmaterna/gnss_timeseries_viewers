# March 1, 2019
# Get the stations within a certain box and write them in the simplest format. 
# For use in easy afterslip calculation. 

import stations_within_radius
import gps_io_functions

coord_box = [-125, -121, 38.5, 42];
infile="../../GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt"
outfile="GPS_ll.txt";

stations = stations_within_radius.get_stations_within_box(coord_box);
[velfield] = gps_io_functions.read_pbo_vel_file(infile);

ofile=open(outfile,'w');
already_added=[];
for i in range(len(velfield.name)):
	if velfield.name[i] in stations:
		if velfield.name[i] in already_added:
			continue;
		else:
			ofile.write("%f %f %s \n" % (velfield.elon[i], velfield.nlat[i], velfield.name[i]) );
			already_added.append(velfield.name[i]);
ofile.close();