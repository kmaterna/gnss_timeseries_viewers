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
Lons
dEv
nd_ll
sig_dNv
dUv
dNv
sig_dEv
"""

import numpy as np 
import scipy.io
import sys

def write_gmt(mat, output_file, output_file_model):
	print("Making text files for %s " % (output_file) );
	print(np.shape(mat['slip_rate']));
	print(np.shape(mat['Kern']));
	print(np.shape(mat['el']));
	print(np.shape(mat['Lons']));

	slip_rate=mat['slip_rate'];
	model_fit = np.dot(mat['Kern'],mat['slip_rate']);
	print(np.shape(model_fit));
	
	# ofile_model = open(output_file_model, 'w');
	# for i in range(len(mat['Lons'])):
	# 	ofile_model.write("%f %f %f %f\n" % (mat['Lons'][i], mat['Lats'][i], model_fit[i*3], model_fit[i*3+1] ) )
	# ofile_model.close();

	ofile=open(output_file,'w');
	for i in range(len(slip_rate)):
		ofile.write(">-Z%f\n" % (slip_rate[i]*1000) );  # in mm
		first_index=mat['el'][i][0]-1;
		second_index=mat['el'][i][1]-1;
		third_index=mat['el'][i][2]-1;
		limit=42;
		if mat['nd_ll'][first_index][1]>limit or mat['nd_ll'][second_index][1]>limit or mat['nd_ll'][third_index][1]>limit:
			continue;
		ofile.write("%f %f\n" % (mat['nd_ll'][first_index][0], mat['nd_ll'][first_index][1]) ) # the first point in the triangle 
		ofile.write("%f %f\n" % (mat['nd_ll'][second_index][0], mat['nd_ll'][second_index][1]) ) # the second point in the triangle
		ofile.write("%f %f\n" % (mat['nd_ll'][third_index][0], mat['nd_ll'][third_index][1]) ) # the third point in the triangle
	ofile.close();
	return;

if __name__=="__main__":

	matlab_file, output_file, output_file_model = '2014_inversion_het_oversmooth.mat', '2014_over.txt', '2014_over_model.txt';
	mat = scipy.io.loadmat(matlab_file);
	write_gmt(mat, output_file, output_file_model);

	matlab_file, output_file, output_file_model = '2014_inversion_het_undersmooth.mat', '2014_under.txt', '2014_under_model.txt';
	mat = scipy.io.loadmat(matlab_file);
	write_gmt(mat, output_file, output_file_model);

	matlab_file, output_file, output_file_model = '2016_inversion_het_oversmooth.mat', '2016_over.txt', '2016_over_model.txt';
	mat = scipy.io.loadmat(matlab_file);
	write_gmt(mat, output_file, output_file_model);

	matlab_file, output_file, output_file_model = '2016_inversion_het_undersmooth.mat', '2016_under.txt', '2016_under_model.txt';
	mat = scipy.io.loadmat(matlab_file);
	write_gmt(mat, output_file, output_file_model);	



