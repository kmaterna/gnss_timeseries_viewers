"""
Make a rotation velocity field

"""

import numpy as np
import haversine


center_of_rotation=[40, -120];  # lat, lon of the center of the rotation. 
omega = 1e-5; # rotation rate, radians per year

input_file="mend_gps_vectors.gmt";
output_file=open("rotation_gps_vectors.gmt",'w')
# NOTE: issues may arise if BM10 and BM1R both exist in the velocity field. 
# Later I will go make sure the input file is clean. 
# Example: I kept P327 and discarded P793, since the uncertainties were lower. 
# 236.426869638 40.4787543283 -6.73 9.03 0.5 0.19 0 0 1 P793_GPS
# 236.426935187 40.4788612703 -7.1 8.85 0.38 0.17 0 0 1 P327_GPS

[lon, lat, vel_east, vel_north, s_east, s_north] = np.loadtxt(input_file,usecols=(0,1,2,3,4,5), unpack=True);

for i in range(len(lon)):
	radius = haversine.distance(center_of_rotation,[lat[i],lon[i]]);
	bearing= haversine.calculate_initial_compass_bearing((center_of_rotation[0],center_of_rotation[1]),(lat[i],lon[i]));
	azimuth=90-bearing;
	radial_velocity = radius*1000*omega;
	new_east = radial_velocity*np.sin(np.deg2rad(azimuth));  # make a rotation velocity field
	new_north = radial_velocity*-1*np.cos(np.deg2rad(azimuth));  # make a rotation velocity field
	output_file.write("%f %f %f %f %f %f 0 0 1 SGPS\n" % (lon[i], lat[i], new_east, new_north, s_east[i], s_north[i]) );



