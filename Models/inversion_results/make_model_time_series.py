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
	folder="HET/nmt_lssq_NA/"
	year="2016"
	matlab_file=folder+year+"_inversion_het.mat";
	output_file=folder+year+"_inversion.txt";
	output_file_model=folder+year+"_model.txt";
	print("Making text files for %s %s" % (folder, year) );
	return matlab_file, output_file, output_file_model;

def read_matlab_file(matlab_file):
	mat = scipy.io.loadmat(matlab_file);
	return mat;

def write_gmt(mat, output_file, output_file_model):

	model_fit = np.dot(mat['Kern'],mat['slip_rate']);
	ofile_model = open(output_file_model, 'w');
	for i in range(len(mat['Lons'])):
		ofile_model.write("%f %f %f %f\n" % (mat['Lons'][i], mat['Lats'][i], model_fit[i*3], model_fit[i*3+1] ) )
	ofile_model.close();

	ofile=open(output_file,'w');
	for i in range(len(mat['slip_rate'])):
		ofile.write(">-Z%f\n" % mat['slip_rate'][i]);
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
	matlab_file, output_file, output_file_model=configure_and_read();
	mat = read_matlab_file(matlab_file);
	write_gmt(mat, output_file, output_file_model);

