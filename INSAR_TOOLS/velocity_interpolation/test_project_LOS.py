# test_project_LOS.py
import do_interpolation
import numpy as np 


def test_project_LOS():

	# Test case 1: 14 degrees ascending orbit, with low-angle camera view, and deformation that's optimally oriented. 
	[Ue, Un, Uu] = [np.cos(np.deg2rad(14)), np.sin(np.deg2rad(14)), 0.0];
	flight_angle = 360 - 14;   # clockwise from north
	incidence_angle = 90;  # horizontal
	expected_answer = 1.0;  # Full projection into the LOS. 
	d_los = do_interpolation.project_to_LOS(Ue, Un, Uu, flight_angle, incidence_angle);
	print("Test case 1: flat incidence angle, def. perpendicular to flight")
	print("computed answer: %f \nexpected answer: %f \n\n" % (d_los, expected_answer) );



	# Test case 2: 14 degrees ascending orbit, with low-angle camera view, and deformation that's parallel to flight (should see nothing). 
	[Ue, Un, Uu] = [-np.sin(np.deg2rad(14)), np.cos(np.deg2rad(14)), 0.0];
	flight_angle = 360 - 14;   # clockwise from north
	incidence_angle = 90;  # horizontal
	expected_answer = 0.0;  # Full projection away from LOS.  
	d_los = do_interpolation.project_to_LOS(Ue, Un, Uu, flight_angle, incidence_angle);
	print("Test case 2: flat incidence angle, def. parallel to flight")
	print("computed answer: %f \nexpected answer: %f \n\n" % (d_los, expected_answer) );


	# Test case 3: 14 degrees ascending orbit, with regular camera view, and deformation 45 degrees to flight (should see small things). 
	[Ue, Un, Uu] = [np.cos(np.deg2rad(59)), np.sin(np.deg2rad(59)), 0.0];
	flight_angle = 360 - 14;   # clockwise from north
	incidence_angle = 30;  # horizontal
	expected_answer = 0.3;  # Full projection into the LOS. 
	d_los = do_interpolation.project_to_LOS(Ue, Un, Uu, flight_angle, incidence_angle);
	print("Test case 3: 30 degree incidence angle, def. 45 degrees to flight")
	print("computed answer: %f \nexpected answer: %f \n\n" % (d_los, expected_answer) );


	# Test case 4: 14 degrees ascending orbit, with regular camera view, and pure vertical deformation (should see large things). 
	[Ue, Un, Uu] = [0 , 0, 1.0];
	flight_angle = 360 - 14;   # clockwise from north
	incidence_angle = 30;  # horizontal
	expected_answer = -0.8;  # Full projection into the LOS. 
	d_los = do_interpolation.project_to_LOS(Ue, Un, Uu, flight_angle, incidence_angle);
	print("Test case 4: 30 degree incidence angle, pure vertical deformation")
	print("computed answer: %f \nexpected answer: %f \n\n" % (d_los, expected_answer) );


	# Test case 5: 14 degrees ascending orbit, with regular camera view, and deformation perpendicular to flight (should see medium things). 
	[Ue, Un, Uu] = [np.cos(np.deg2rad(14)), np.sin(np.deg2rad(14)), 0.0];
	flight_angle = 360 - 14;   # clockwise from north
	incidence_angle = 30;  # horizontal
	expected_answer = 0.5;  # Full horizontal projection into the LOS. 
	d_los = do_interpolation.project_to_LOS(Ue, Un, Uu, flight_angle, incidence_angle);
	print("Test case 5: 30 degree incidence angle, def. 45 degrees to flight")
	print("computed answer: %f \nexpected answer: %f \n\n" % (d_los, expected_answer) );


	return;

if __name__=="__main__":
	test_project_LOS();

