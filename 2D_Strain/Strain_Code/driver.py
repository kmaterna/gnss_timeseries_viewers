"""
Driver program for strain calculation

"""

import common_io_functions
import delaunay_strain
import hammond_strain
import gpsgridder_strain
import visr_strain


def driver_1d(strain_method):
	[MyParams] = common_io_functions.configure(strain_method);
	[myVelfield] = common_io_functions.inputs(MyParams);
	[xdata, ydata, polygon_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11] = compute_dict[strain_method](myVelfield, MyParams);
	common_io_functions.outputs_1d(xdata, ydata, polygon_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, myVelfield, MyParams);	
	return;

def driver_2d(strain_method):
	[MyParams] = common_io_functions.configure(strain_method);
	[myVelfield] = common_io_functions.inputs(MyParams);
	[xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11] = compute_dict[strain_method](myVelfield, MyParams);
	common_io_functions.outputs_2d(xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, myVelfield, MyParams);	
	return;



# Are we producing gridded interpolated data, or static polygons?
driver_dict={
	"delaunay":driver_1d, 
	"hammond":driver_1d,
	"gpsgridder":driver_2d,
	"visr": driver_2d };

# Where does the compute method live? 
compute_dict={
	"delaunay":delaunay_strain.compute, 
	"hammond":hammond_strain.compute,
	"gpsgridder":gpsgridder_strain.compute,
	"visr":visr_strain.compute };



if __name__=="__main__":

	strain_method="gpsgridder"
	driver_dict[strain_method](strain_method);

	strain_method="delaunay"
	driver_dict[strain_method](strain_method);

	strain_method="hammond"
	driver_dict[strain_method](strain_method);

	strain_method="visr"
	driver_dict[strain_method](strain_method);

