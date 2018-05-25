# The configure, input, and output steps for GPS Strain analysis. 
# map_range=[-125, -121, 37.0, 42.2]; # Northern California
# map_range=[-121, -115, 32, 36]; # Southern California
# map_range=[-123, -118, 34, 38]; # Central California
# map_range=[-125, -114, 32.0, 42.0]; # ALL California	

import numpy as np 
import subprocess, sys, os
import collections
import gps_io_functions
import netcdf_io_functions
import delaunay_strain
import gpsgridder_strain
import hammond_strain



Params=collections.namedtuple("Params",['strain_method','input_file','map_range','coord_box','num_years','max_sigma','grid_inc','outdir','gmtfile']);



# ----------------- CONFIGURE -------------------------
def configure(strain_method):
	input_file="../../GPS_POS_DATA/PBO_Velocity_Files/NAM08_pbovelfile_feb2018.txt";
	map_range=[-125, -121, 37.0, 42.2]; # Northern California
	map_range_string = str(map_range[0])+'/'+str(map_range[1])+'/'+str(map_range[2])+'/'+str(map_range[3]);
	num_years=3.0;
	max_sigma=2.0;
	[grid_inc, coord_box, outdir, gmtfile] = get_tunable_options(strain_method, map_range);
	MyParams=Params(strain_method=strain_method, input_file=input_file, map_range=map_range_string, coord_box=coord_box, num_years=num_years, max_sigma=max_sigma, grid_inc=grid_inc, outdir=outdir, gmtfile=gmtfile);
	return [MyParams];


# A couple of options that change based on strain method. 
def get_tunable_options(strain_method, map_range):
	if strain_method=="gpsgridder":
		grid_inc   =0.02;
		coord_box  =[map_range[0]-1, map_range[1]+3, map_range[2]-2, map_range[3]+2];
		outdir     ="../GPSgridder/";
		gmtfile    ="gpsgridder_gmt.gmt";

	elif strain_method=="delaunay":
		grid_inc   =0.04; # larger interval for convenience, because it's slow. 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		outdir     ="../Delaunay/"
		gmtfile    ="delaunay_gmt.gmt"

	elif strain_method=="hammond":
		grid_inc   =0.04; # larger interval for convenience, because it's slow. 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		outdir     ="../Hammond/"
		gmtfile    ="hammond_gmt.gmt"

	else:
		print("ERROR: "+strain_method+" is not a known strain method. ");
		sys.exit(1);

	return [grid_inc, coord_box, outdir, gmtfile];



# ----------------- INPUTS -------------------------
def inputs(MyParams):
	# Purpose: generate input velocity field. 
	[myVelfield]=gps_io_functions.read_pbo_vel_file(MyParams.input_file);  # read the raw velfield from file. 
	print(len(myVelfield.name));
	[myVelfield]=gps_io_functions.clean_velfield(myVelfield, MyParams.num_years, MyParams.max_sigma, MyParams.coord_box);
	print(len(myVelfield.name));
	[myVelfield]=gps_io_functions.remove_duplicates(myVelfield);
	print(len(myVelfield.name));		
	return [myVelfield];





# ----------------- OUTPUTS -------------------------

def outputs_2d(xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, myVelfield, MyParams):
	print("Sending outputs")
	outfile=open(MyParams.outdir+"tempgps.txt",'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0.0\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i]) );
	outfile.close();
	netcdf_io_functions.produce_output_netcdf(xdata, ydata, I2nd, 'per yr', MyParams.outdir+'I2nd.nc');
	netcdf_io_functions.flip_if_necessary(MyParams.outdir+'I2nd.nc');
	netcdf_io_functions.produce_output_netcdf(xdata, ydata, rot, 'per yr', MyParams.outdir+'rot.nc');
	netcdf_io_functions.flip_if_necessary(MyParams.outdir+'rot.nc');	
	write_grid_eigenvectors(xdata, ydata, e1, e2, v00, v01, v10, v11, MyParams);
	os.chdir(MyParams.outdir)
	subprocess.call("../Strain_Code/"+MyParams.gmtfile+" "+MyParams.map_range,shell=True);
	return;

def write_grid_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11, MyParams):
	# Need eigs_interval and outdir from MyParams. 
	positive_file=open(MyParams.outdir+"positive_eigs.txt",'w');
	negative_file=open(MyParams.outdir+"negative_eigs.txt",'w');
	eigs_dec=12;

	for j in range(len(ydata)):
		for k in range(len(xdata)):
			if np.mod(j,eigs_dec)==0 and np.mod(k,eigs_dec)==0:
				if w1[j][k]>0:
					scale=w1[j][k];
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k]*scale, v10[j][k]*scale) );
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k]*scale, -v10[j][k]*scale) );
				if w1[j][k]<0:
					scale=w1[j][k];
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k]*scale, v10[j][k]*scale) );
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k]*scale, -v10[j][k]*scale) );
				if w2[j][k]>0:
					scale=w2[j][k];
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k]*scale, v11[j][k]*scale) );
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k]*scale, -v11[j][k]*scale) );
				if w2[j][k]<0:
					scale=w2[j][k];
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k]*scale, v11[j][k]*scale) );
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k]*scale, -v11[j][k]*scale) );
	return;






def outputs_1d(xcentroid, ycentroid, polygon_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, myVelfield, MyParams):

	rotfile=open(MyParams.outdir+"rotation.txt",'w');
	I2ndfile=open(MyParams.outdir+"I2nd.txt",'w');
	positive_file=open(MyParams.outdir+"positive_eigs.txt",'w');
	negative_file=open(MyParams.outdir+"negative_eigs.txt",'w');

	outfile=open(MyParams.outdir+"tempgps.txt",'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0.0\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i]) );
	outfile.close();

	for i in range(len(I2nd)):
		# Write the triangle's rotation
		rotfile.write("> -Z"+str(rot[i])+"\n");
		rotfile.write(str(polygon_vertices[i,0,0])+" "+str(polygon_vertices[i,0,1])+"\n");
		rotfile.write(str(polygon_vertices[i,1,0])+" "+str(polygon_vertices[i,1,1])+"\n");
		rotfile.write(str(polygon_vertices[i,2,0])+" "+str(polygon_vertices[i,2,1])+"\n");

		# Write the triangle's I2
		I2ndfile.write("> -Z"+str(I2nd[i])+"\n");
		I2ndfile.write(str(polygon_vertices[i,0,0])+" "+str(polygon_vertices[i,0,1])+"\n");
		I2ndfile.write(str(polygon_vertices[i,1,0])+" "+str(polygon_vertices[i,1,1])+"\n");
		I2ndfile.write(str(polygon_vertices[i,2,0])+" "+str(polygon_vertices[i,2,1])+"\n");

		# Write the eigenvectors and eigenvalues
		write_single_eigenvector(positive_file, negative_file, e1[i], v00[i], v10[i], xcentroid[i], ycentroid[i]);
		write_single_eigenvector(positive_file, negative_file, e2[i], v01[i], v11[i], xcentroid[i], ycentroid[i]);
	
	print(max(I2nd));
	print(max(rot));
	print(min(rot));

	rotfile.close();
	I2ndfile.close();
	positive_file.close();
	negative_file.close();
	os.chdir(MyParams.outdir)	
	subprocess.call("../Strain_Code/"+MyParams.gmtfile+" "+MyParams.map_range,shell=True);

	return;


def write_single_eigenvector(positive_file, negative_file, e, v0, v1, x, y):
	# e = eigenvalue, [v0, v1] = eigenvector. 
	# Writes a single eigenvector eigenvalue pair. 
	# Also has functionality to saturate eigenvectors so they don't blow up. 
	overall_max=40.0;
	scale=0.4*e;

	vx=v0*scale;
	vy=v1*scale;
	if np.sqrt(vx*vx+vy*vy)>overall_max:
		scale=scale*(overall_max/np.sqrt(vx*vx+vy*vy))
		vx=v0*scale;
		vy=v1*scale;

	if e>0:
		positive_file.write("%s %s %s %s 0 0 0\n" % (x, y, vx, vy ));
		positive_file.write("%s %s %s %s 0 0 0\n" % (x, y, -vx, -vy ));
	else:	
		negative_file.write("%s %s %s %s 0 0 0\n" % (x, y, vx, vy ));
		negative_file.write("%s %s %s %s 0 0 0\n" % (x, y, -vx, -vy ));
	return;



