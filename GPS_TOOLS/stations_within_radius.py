# Python functions to take an input coordinate, and 
# determine which PBO stations are within a certain radius of that point. 


import haversine
import gps_io_functions

# Reference: Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);

def get_stations_within_radius(center, radius, coord_box=[]):
	[input_file, center, radius, num_years, max_sigma, coord_box] = configure(center, radius, coord_box);
	myVelfield = inputs(input_file, num_years, max_sigma, coord_box);
	close_stations, rad_distance = compute(myVelfield, center, radius);	
	return close_stations, rad_distance;

def configure(center, radius, coord_box):
	input_file="../GPS_POS_DATA/PBO_Velocity_Files/NAM08_pbovelfile_feb2018.txt";
	num_years=3.0;
	max_sigma=2.0;
	if coord_box==[]:
		coord_box=[-126, -120, 36.0, 43]; # Northern California
	return [input_file, center, radius, num_years, max_sigma, coord_box];


def inputs(input_file, num_years, max_sigma, coord_box):
	# Purpose: generate input velocity field. 
	[myVelfield]=gps_io_functions.read_pbo_vel_file(input_file);  # read the raw velfield from file. 
	[myVelfield]=gps_io_functions.clean_velfield(myVelfield, num_years, max_sigma, coord_box);
	[myVelfield]=gps_io_functions.remove_duplicates(myVelfield);
	return myVelfield;


def compute(myVelfield, center, radius):
	close_stations=[];
	rad_distance=[];
	for i in range(len(myVelfield.name)):
		mydist = haversine.distance([center[1], center[0]],[myVelfield.nlat[i],myVelfield.elon[i]])
		if mydist <= radius:
			rad_distance.append(mydist);
			close_stations.append(myVelfield.name[i]);
	print("Returning %d stations that are within %f km of center %f, %f" % (len(close_stations), radius, center[0], center[1]) )
	return close_stations, rad_distance;


if __name__=="__main__":
	center=[-124,40];
	radius=120;
	close_stations, rad_distance = get_stations_within_radius(center, radius);
