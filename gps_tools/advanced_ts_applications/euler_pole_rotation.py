# Take a velocity (or list of velocities) in one reference frame, 
# And rotate them into a second reference frame using the Euler Pole of 
# that reference frame transformation.
import gps_tools.file_io.io_nota
import gps_tools.file_io.io_other
import numpy as np 
import matplotlib.pyplot as plt 
import collections
import gps_io_functions
Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);

def configure():
	input_file="../GPS_POS_DATA/Velocity_Files/IGS08_pbovelfile_feb2018.txt";
	output_file="../GPS_POS_DATA/testing_euler_pole/nam.txt"
	# North America ITRF97 Euler Pole (https://web.ics.purdue.edu/~ecalais/teaching/kinematics/plate_motions.pdf)
	EP_lon=-79.08;   # degrees east
	EP_lat=-2.39;   # degrees north
	EP_rot=0.199;  # degrees per million years
	euler_pole=[EP_lon, EP_lat, EP_rot];  
	return [input_file, output_file, euler_pole];

def configure_SNGV():
	input_file="../GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt";
	output_file="../GPS_POS_DATA/testing_euler_pole/SNGV.txt"
	# North America SNGV euler pole (https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/1998TC001088).
	# Dixon et al. 2000
	EP_lon=-137;   # degrees east
	EP_lat=17;   # degrees north
	EP_rot=0.28;  # degrees per million years
	euler_pole=[EP_lon, EP_lat, EP_rot];  
	return [input_file, output_file, euler_pole];


def compute_euler_pole_rotation(vels, euler_pole):
	# This is where the important math lives. 
	# Operates on a velocity object, and an euler pole (3-item list)
	# Returns another velocity object
	new_east=[]; new_north=[]; new_up=[];

	for i in range(len(vels.name)):

		R_point = get_r(vels.elon[i], vels.nlat[i]);
		R_ep = get_r(euler_pole[0], euler_pole[1]);
		unit_ep = get_unit_vector(R_ep);
		omega_raw = degma2radyr(euler_pole[2]);
		omega = omega_raw * unit_ep; # in radians per year

		velocity_of_transformation = np.cross(omega, R_point);  # velocity at the station from the euler pole rotation
		velocity_of_transformation = velocity_of_transformation*1000;  # mm/yr in x, y, z
		
		xvel = velocity_of_transformation[0];
		yvel = velocity_of_transformation[1];
		zvel = velocity_of_transformation[2];
		[east_transform, north_transform, up_transform] = xyz2enu(xvel, yvel, zvel, vels.elon[i], vels.nlat[i]);

		east_after = vels.e[i]-east_transform;
		north_after = vels.n[i]-north_transform;
		up_after = vels.u[i]-up_transform; 
		new_east.append(east_after);
		new_north.append(north_after);
		new_up.append(up_after);

	new_vels=Velfield(name=vels.name, nlat=vels.nlat, elon=vels.elon, n=new_north, e=new_east, u=new_up, sn=vels.sn, se=vels.sn, su=vels.su, first_epoch=vels.first_epoch, last_epoch=vels.last_epoch);
	return new_vels;


def get_unit_vector(vec):
	mag=np.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	vec=np.divide(vec,mag);
	return vec;

def degma2radyr(omega):
	# degrees/Ma to radians/yr
	radyr = omega*(np.pi/180)*1e-6;
	return radyr;

def get_r(lon,lat):
	# Vector from center of earth to the point in question
	# Definitions: 
	# 0 longitude: x = 0
	# 0 latitude: z = 0
	# North is positive z. 
	R_fixed = 6378000;  # In meters
	R_equatorial_disk = R_fixed * np.cos(np.deg2rad(lat)); 
	T_equatorial_disk = np.deg2rad(lon); 
	X = R_equatorial_disk*np.cos(T_equatorial_disk); 
	Y = R_equatorial_disk*np.sin(T_equatorial_disk); 
	Z = np.sqrt(R_fixed*R_fixed - X*X - Y*Y); 
	if lat<0:
		Z=Z*-1;
	return [X, Y, Z];

def get_unit_east(lon, lat):
	T_equatorial_disk = np.deg2rad(lon); 
	x = -np.sin(T_equatorial_disk);
	y = np.cos(T_equatorial_disk);
	return [x, y, 0];

def xyz2enu(x, y, z, lon, lat):
	vel_vector = [x, y, z];
	# Convert velocities from xyz to east north up
	# Assuming spherical earth
	# Assuming horizontal movement only 
	# The unit east vector is the only one we need. 
	# Then we take the dot product with the unit east vector. 
	# The north component is the remainder. 
	unit_east = get_unit_east(lon, lat);
	e = np.dot(vel_vector,unit_east);
	n = np.sqrt(x*x + y*y + z*z - e*e);
	if z < 0:
		n=n*-1;
	u = 0;
	return [e, n, u];


def plot_na_itrf(origvel, aftervel):
	plt.quiver(origvel.elon, origvel.nlat, origvel.e, origvel.n, scale=500,color='black');
	plt.quiver(aftervel.elon, aftervel.nlat, aftervel.e, aftervel.n, scale=500,color='red');
	plt.xlim([-125, -120]);
	plt.ylim([38, 42]);
	plt.legend(['NA','SNGV'])
	plt.savefig('../GPS_POS_DATA/testing_euler_pole/rotations.eps');
	return;


if __name__=="__main__":
	# A test situation, converting ITRF into NA
	# [input_file, output_file, euler_pole] = configure();
	[input_file, output_file, euler_pole] = configure_SNGV();
	[vels] = gps_tools.file_io.io_nota.read_pbo_vel_file(input_file);
	new_vels = compute_euler_pole_rotation(vels, euler_pole);
	# gps_io_functions.write_humanread_vel_file(new_vels, output_file);
	gps_tools.file_io.io_other.write_gmt_velfile(new_vels, output_file);
	plot_na_itrf(vels, new_vels);

