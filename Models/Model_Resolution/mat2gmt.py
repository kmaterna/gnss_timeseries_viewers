# April 2019
# Read the matlab file
# Write a text file
# Plot in GMT

"""
stations_used
Lats
slip_rate
el
nd
Kern
sig_dUv
__version__
Lons
__globals__
d
dEv
nd_ll
sig_dNv
dUv
dNv
sig_dEv
"""

import numpy as np 
import scipy.io

def configure_and_read():
	matlab_file="resolution_for_Kathryn.mat";
	ll_file="2014_inversion_het.mat"
	output_file="resolution.txt";
	print("Making text files for %s " % (matlab_file) );
	return matlab_file, ll_file, output_file;

def read_matlab_file(matlab_file):
	mat = scipy.io.loadmat(matlab_file);
	return mat;

def write_gmt(mat, ll, output_file):

	# model_fit = np.dot(mat['Kern'],mat['slip_rate']);
	# ofile_model = open(output_file_model, 'w');
	# for i in range(len(mat['Lons'])):
	# 	ofile_model.write("%f %f %f %f\n" % (mat['Lons'][i], mat['Lats'][i], model_fit[i*3], model_fit[i*3+1] ) )
	# ofile_model.close();

	print(np.shape(ll['slip_rate']));
	print(np.shape(mat['mean_Kern']));
	mean_Kern=np.array(mat['mean_Kern']);
	mean_Kern=mean_Kern.T;

	ofile=open(output_file,'w');
	for i in range(len(mean_Kern)):
		ofile.write(">-Z%f\n" % (mean_Kern[i]*1000) );  # in mm
		first_index=ll['el'][i][0]-1;
		second_index=ll['el'][i][1]-1;
		third_index=ll['el'][i][2]-1;
		limit=42;
		if ll['nd_ll'][first_index][1]>limit or ll['nd_ll'][second_index][1]>limit or ll['nd_ll'][third_index][1]>limit:
			continue;
		ofile.write("%f %f\n" % (ll['nd_ll'][first_index][0], ll['nd_ll'][first_index][1]) ) # the first point in the triangle 
		ofile.write("%f %f\n" % (ll['nd_ll'][second_index][0], ll['nd_ll'][second_index][1]) ) # the second point in the triangle
		ofile.write("%f %f\n" % (ll['nd_ll'][third_index][0], ll['nd_ll'][third_index][1]) ) # the third point in the triangle
	ofile.close();
	return;

if __name__=="__main__":
	matlab_file, ll_file, output_file=configure_and_read();
	mat = read_matlab_file(matlab_file);
	ll = read_matlab_file(ll_file);
	write_gmt(mat, ll, output_file);

