# Take a velocity (or list of velocities) in one reference frame, 
# And rotate them into a second reference frame using the Euler Pole of 
# that reference frame transformation. 
import gps_io_functions
import haversine

def configure():
	input_file="../GPS_POS_DATA/Velocity_Files/IGS08_pbovelfile_feb2018.txt";
	output_file="../GPS_POS_DATA/testing_euler_pole/nam.txt"
	# North America ITRF97 Euler Pole (https://web.ics.purdue.edu/~ecalais/teaching/kinematics/plate_motions.pdf)
	EP_lon=-79.08;   # degrees east
	EP_lat=-2.39;   # degrees north
	EP_rot=0.199;  # degrees per million years
	euler_pole=[EP_lon, EP_lat, EP_rot];  
	return [input_file, output_file, euler_pole];

def compute_euler_pole_rotation(vels, euler_pole):
	# This is where the important math lives. 

	
	return vels;


if __name__=="__main__":
	# A test situation, converting ITRF into NA
	[input_file, output_file, euler_pole] = configure();
	[vels] = gps_io_functions.read_pbo_vel_file(input_file);
	new_vels = compute_euler_pole_rotation(vels, euler_pole);
	gps_io_functions.write_humanread_vel_file(new_vels, output_file);