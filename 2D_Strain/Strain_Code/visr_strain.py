# Use the interpolation scheme of Zheng-Kang Shen et al.
# Shen, Z.-K., M. Wang, Y. Zeng, and F. Wang, Strain determination using spatially discrete geodetic data, 
# Bull. Seismol. Soc. Am., 105(4), 2117-2127, doi: 10.1785/0120140247, 2015.
# http://scec.ess.ucla.edu/~zshen/visr/visr.html

# The fortran files must be compiled and linked like this: 
# gfortran -c voronoi_area_version.f90
# gfortran visr.f voronoi_area_version.o -o visr.exe


import numpy as np 
import subprocess
import strain_tensor_toolbox


def compute(myVelfield, MyParams):

	# write_into_fortran_format();
	# call_fortran_compute();

	# Here we will make the strain compute function, using visr's fortran code. 
	# It will output a large text file. 
	# We convert that text file into grids, which we will write as GMT grd files. 

	[xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11] = make_output_grids_from_strain_out("visr/strain.out");
	return [xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11];


def make_output_grids_from_strain_out(infile):
	ifile=open(infile,'r');
	x=[]; y=[]; rotation=[]; exx=[]; exy=[]; eyy=[];
	for line in ifile:
		temp=line.split();
		if 'index' in temp or 'longitude' in temp or 'deg' in temp:
			continue;
		else:
			x.append(float(temp[0]));
			y.append(float(temp[1]));
			rotation.append(float(temp[7]));
			exx.append(float(temp[9]));
			exy.append(float(temp[11]));
			eyy.append(float(temp[13]));
	ax1=set(x);
	ax2=set(y);
	xlen=len(ax1);
	ylen=len(ax2);

	xaxis=sorted(ax1);
	yaxis=sorted(ax2);
	
	# Loop through x and y lists, find the index of the coordinates in the xaxis and yaxis sets, 
	# Then place them into the 2d arrays. 
	# Then go compute I2nd, eigenvectors and eigenvalues. 

	I2nd=np.zeros((ylen, xlen)); max_shear=np.zeros((ylen, xlen)); rot=np.zeros((ylen, xlen)); e1=np.zeros((ylen, xlen)); e2=np.zeros((ylen, xlen));
	v00=np.zeros((ylen, xlen)); v01=np.zeros((ylen, xlen)); v10=np.zeros((ylen, xlen)); v11=np.zeros((ylen, xlen)); 
	print(np.shape(xaxis));
	print(np.shape(I2nd))

	for i in range(len(x)):
		xindex=xaxis.index(x[i]);
		yindex=yaxis.index(y[i]);
		rot[yindex][xindex]=rotation[i];
		I2nd[yindex][xindex] = np.log10(np.abs(strain_tensor_toolbox.second_invariant(exx[i], exy[i], eyy[i])));
		[e11, e22, v1] = strain_tensor_toolbox.eigenvector_eigenvalue(exx[i], exy[i], eyy[i]);
		e1[yindex][xindex]= -e11;
		e2[yindex][xindex]= -e22;
		v00[yindex][xindex]=v1[0][0];
		v10[yindex][xindex]=v1[0][1];
		v01[yindex][xindex]=v1[1][0];
		v11[yindex][xindex]=v1[1][1];
		max_shear[yindex][xindex] = (e11 - e22)/2;


	return [xaxis, yaxis, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11];


