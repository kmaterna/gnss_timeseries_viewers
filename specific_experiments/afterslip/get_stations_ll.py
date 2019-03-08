# March 1, 2019
# Get the stations within a certain box and write them in the simplest format. 
# For use in easy afterslip calculation. 

import stations_within_radius
import gps_io_functions
import gps_input_pipeline
import datetime as dt 

coord_box = [-125, -120, 42, 46.0];
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

			ofile.write("%f %f %s \n" % (velfield.elon[i], velfield.nlat[i], velfield.name[i]) );  # format for afterslip
			
			# # IF doing this for GRACE: 
			# [newData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(velfield.name[i], 'pbo');			
			# start = float(dt.datetime.strftime(newData.dtarray[0],"%Y"))-1.0;
			# end = float(dt.datetime.strftime(newData.dtarray[-1],"%Y"))+1.0;
			# if start <= 2002: 
			# 	start=2002; 
			# if end >= 2018:
			# 	end=2018;
			# ofile.write("%s %f %f %.2f %.2f \n" % (newData.name, newData.coords[0], newData.coords[1], start, end ) );  # format for grace
			
			already_added.append(velfield.name[i]);
ofile.close();