"""
Take a set of velocities, establish delaunay triangles, 
solve a linear inversion problem for the components of the velocity gradient tensor
at the centroid of each triangle. 
The strain rate tensor and the rotation tensor can be readily computed 
from the symmetric and anti-symmetric parts of the velocity gradient tensor. 
Plot the outputs. 

Following a technique learned in Brad Hagar's geodynamics class, and 
modeled off of advice from 2007 Journal of Geodynamcis paper:
ftp://ftp.ingv.it/pub/salvatore.barba/RevEu/Cai_StrainBIFROST_2007.pdf
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy.spatial import Delaunay
from numpy.linalg import inv, det, eig
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def configure():
	input_file="mend_gps_vectors.gmt";
	# input_file="rotation_gps_vectors.gmt";
	# NOTE: issues may arise if BM10 and BM1R both exist in the velocity field. 
	# Later I will go make sure the input file is clean. 
	# Example: I kept P327 and discarded P793, since the uncertainties were lower. 
	# 236.426869638 40.4787543283 -6.73 9.03 0.5 0.19 0 0 1 P793_GPS
	# 236.426935187 40.4788612703 -7.1 8.85 0.38 0.17 0 0 1 P327_GPS
	return [input_file];


def inputs(input_file):
	[lon, lat, vel_east, vel_north, s_east, s_north] = np.loadtxt(input_file,usecols=(0,1,2,3,4,5), unpack=True);
	return [lon-360, lat, vel_east, vel_north, s_east, s_north];




def process(lon, lat, vel_east, vel_north, s_east, s_north):
	z = np.array([lon,lat]);
	z = z.T;
	tri=Delaunay(z);

	triangle_vertices = z[tri.simplices];
	trishape = np.shape(triangle_vertices);  # 516 x 3 x 2, for example
	print trishape[0]

	# We are going to solve for the velocity gradient tensor at the centroid of each triangle. 
	centroids=[];
	for i in range(trishape[0]):
		xcor_mean = np.mean([triangle_vertices[i,0,0],triangle_vertices[i,1,0],triangle_vertices[i,2,0]]);
		ycor_mean = np.mean([triangle_vertices[i,0,1],triangle_vertices[i,1,1],triangle_vertices[i,2,1]]);
		centroids.append([xcor_mean,ycor_mean]);
	xcentroid=[x[0] for x in centroids];
	ycentroid=[x[1] for x in centroids];

	dVEdE=np.zeros(trishape[0]);
	dVEdN=np.zeros(trishape[0]);
	dVNdE=np.zeros(trishape[0]);
	dVNdN=np.zeros(trishape[0]);
	Rdet=np.zeros(trishape[0]);
	Tdet=np.zeros(trishape[0]);
	Varray=[];
	Warray=[];

	# for each triangle:
	for i in range(trishape[0]):

		# Get the velocities of each vertex (VE1, VN1, VE2, VN2, VE3, VN3)
		# Get velocities for Vertex 1 (triangle_vertices[i,0,0] and triangle_vertices[i,0,1])
		xindex1 = np.where(lon==triangle_vertices[i,0,0])
		yindex1 = np.where(lat==triangle_vertices[i,0,1])
		index1=np.intersect1d(xindex1,yindex1);
		xindex2 = np.where(lon==triangle_vertices[i,1,0])
		yindex2 = np.where(lat==triangle_vertices[i,1,1])
		index2=np.intersect1d(xindex2,yindex2);
		xindex3 = np.where(lon==triangle_vertices[i,2,0])
		yindex3 = np.where(lat==triangle_vertices[i,2,1])
		index3=np.intersect1d(xindex3,yindex3);

		VE1=vel_east[index1]; VN1=vel_north[index1];
		VE2=vel_east[index2]; VN2=vel_north[index2];
		VE3=vel_east[index3]; VN3=vel_north[index3];
		obs_vel = np.array([[VE1[0]],[VN1[0]],[VE2[0]],[VN2[0]],[VE3[0]],[VN3[0]]]);

		# Get the distance between centroid and vertex
		dE1 = triangle_vertices[i,0,0]-xcentroid[i];
		dE2 = triangle_vertices[i,1,0]-xcentroid[i];
		dE3 = triangle_vertices[i,2,0]-xcentroid[i];
		dN1 = triangle_vertices[i,0,1]-ycentroid[i];
		dN2 = triangle_vertices[i,1,1]-ycentroid[i];
		dN3 = triangle_vertices[i,2,1]-ycentroid[i];


		Design_Matrix = np.array([[1,0,dE1,dN1,0,0],[0,1,0,0,dE1,dN1],[1,0,dE2,dN2,0,0],[0,1,0,0,dE2,dN2],[1,0,dE3,dN3,0,0],[0,1,0,0,dE3,dN3]]);

		# Invert to get the components of the velocity gradient tensor. 
		DMinv = inv(Design_Matrix);
		vel_grad = np.dot(DMinv, obs_vel);  # this is the money step. 
		VE_centroid=vel_grad[0][0];
		VN_centroid=vel_grad[1][0];
		dVEdE[i]=vel_grad[2][0];
		dVEdN[i]=vel_grad[3][0];
		dVNdE[i]=vel_grad[4][0];
		dVNdN[i]=vel_grad[5][0];

		# The strain rate tensor (not saved for each triangle)
		T = np.array([[dVEdE[i], 0.5*(dVEdN[i]+dVNdE[i])],[0.5*(dVNdE[i]+dVEdN[i]), dVNdN[i]]])

		# The rotation rate tensor (not saved for each triangle)
		R = np.array([[0, 0.5*(dVEdN[i]-dVNdE[i])],[0.5*(dVNdE[i]-dVEdN[i]),0]])

		w,v=np.linalg.eig(T);  # The eigenvectors and eigenvalues (principal strains) of the strain rate tensor
		Warray.append(w);  # eigenvalues
		Varray.append(v);  # eigenvectors

		# The determinant of the rotation rate tensor (for plotting)
		Rdet[i]=np.linalg.det(R);

		# The determinant of the strain rate tensor (the 2nd invariant for the 2D case)
		Tdet[i]=np.linalg.det(T);

	return [z, tri, centroids, Tdet, Rdet, Warray, Varray];




def outputs(lon, lat, vel_east, vel_north, s_east, s_north, z, tri, centroids, Tdet, Rdet,Warray,Varray):
	[lon_ca, lat_ca] = np.loadtxt("california_bdr",unpack=True);
	[lon_nv, lat_nv] = np.loadtxt("nevada_bdr",unpack=True);
	[lon_or, lat_or] = np.loadtxt("oregon_bdr",unpack=True);

	# Cleaning up values that blew up. 
	max_allowable=30;
	for i in range(len(Tdet)):
		if Tdet[i]>max_allowable:
			print Tdet[i];
			Tdet[i]=max_allowable;
			# Warray[i]=[0,0];
	min_allowable=-2000;
	for i in range(len(Tdet)):
		if Tdet[i]<min_allowable:
			print Tdet[i];
			Tdet[i]=min_allowable;
			Warray[i]=[0,0];
	max_allowable=2000;
	for i in range(len(Rdet)):
		if Rdet[i]>max_allowable:
			Rdet[i]=max_allowable;


	# Plot for Tdet
	plt.figure();
	plt.plot(lon,lat,'.k')
	plt.plot(lon_ca,lat_ca,'k',linewidth=3)
	plt.plot(lon_nv,lat_nv,'k')
	plt.plot(lon_or,lat_or,'k')
	xcentroid=[x[0] for x in centroids];
	ycentroid=[x[1] for x in centroids];
	plt.quiver(lon,lat,vel_east,vel_north)
	plt.triplot(z[:,0], z[:,1], tri.simplices.copy())

	# Start to plot the patches
	triangle_vertices = z[tri.simplices];
	trishape = np.shape(triangle_vertices);  # 516 x 3 x 2, for example
	patches=[];
	for i in range(trishape[0]):
		v=np.array([[triangle_vertices[i,0,0], triangle_vertices[i,0,1]],[triangle_vertices[i,1,0], triangle_vertices[i,1,1]],[triangle_vertices[i,2,0], triangle_vertices[i,2,1]]]);
		patches.append(Polygon(v,True));
	p = PatchCollection(patches, alpha=0.4)

	p.set_array(np.array(Tdet))
	ax=plt.gca();
	ax.add_collection(p)
	plt.colorbar(p, ax=ax)
	plt.xlim([-125, -118])
	plt.ylim([36.5, 42.5])
	plt.title("Determinant of Strain Rate Tensor")
	plt.savefig("Tdet.jpg");

	scaling=800.0;
	for i in range(len(Warray)):
		eigvals=Warray[i];
		eigvecs=Varray[i];  # each column of this matrix is an eigenvector
		if eigvals[0]<0:
			print i;
			print eigvals[0]
			print eigvals[1];
			plt.quiver(xcentroid[i],ycentroid[i],eigvals[0]*eigvecs[0][0],eigvals[0]*eigvecs[1][0],color='red',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
			plt.quiver(xcentroid[i],ycentroid[i],-eigvals[0]*eigvecs[0][0],-eigvals[0]*eigvecs[1][0],color='red',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
		if eigvals[1]<0:
			plt.quiver(xcentroid[i],ycentroid[i],eigvals[1]*eigvecs[0][1],eigvals[1]*eigvecs[1][1],color='red',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
			plt.quiver(xcentroid[i],ycentroid[i],-eigvals[1]*eigvecs[0][1],-eigvals[1]*eigvecs[1][1],color='red',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
		if eigvals[0]>0:
			plt.quiver(xcentroid[i],ycentroid[i],eigvals[0]*eigvecs[0][0],eigvals[0]*eigvecs[1][0],color='blue',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
			plt.quiver(xcentroid[i],ycentroid[i],-eigvals[0]*eigvecs[0][0],-eigvals[0]*eigvecs[1][0],color='blue',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
		if eigvals[1]>0:
			plt.quiver(xcentroid[i],ycentroid[i],eigvals[1]*eigvecs[0][1],eigvals[1]*eigvecs[1][1],color='blue',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
			plt.quiver(xcentroid[i],ycentroid[i],-eigvals[1]*eigvecs[0][1],-eigvals[1]*eigvecs[1][1],color='blue',scale=scaling,headaxislength=0,headwidth=1, headlength=0);
	h1,=plt.plot([0,1],[0,1],'r',linewidth=3,label='Contraction');
	h2,=plt.plot([0,1],[0,1],'b',linewidth=3,label='Extension');
	plt.legend(handles=[h1,h2])

	plt.title('Principal Strain Axes in Northern California from GPS')
	plt.xlim([-124.9, -122.6]);
	plt.ylim([39.4, 41.4]);
	plt.savefig("Tdet_closeup.jpg");

	plt.title('Principal Strain Axes in CSZ')
	plt.xlim([-125.3, -122.6]);
	plt.ylim([39.9, 42.4]);
	plt.savefig("Tdet_CSZ.jpg");


	plt.title('Principal Strain Axes in Bay Area')
	plt.xlim([-124.1, -121.6]);
	plt.ylim([37.0, 39]);
	plt.savefig("Tdet_closeup_Bay_Area.jpg");
	plt.close();


	# Plot for Rdet
	plt.figure();
	plt.plot(lon,lat,'.k')
	plt.plot(lon_ca,lat_ca,'k',linewidth=3)
	plt.plot(lon_nv,lat_nv,'k')
	plt.plot(lon_or,lat_or,'k')
	xcentroid=[x[0] for x in centroids];
	ycentroid=[x[1] for x in centroids];
	# plt.plot(xcentroid,ycentroid,'.r')
	plt.quiver(lon,lat,vel_east,vel_north)
	plt.triplot(z[:,0], z[:,1], tri.simplices.copy())

	# Start to plot the patches
	triangle_vertices = z[tri.simplices];
	trishape = np.shape(triangle_vertices);  # 516 x 3 x 2, for example
	patches=[];
	for i in range(trishape[0]):
		v=np.array([[triangle_vertices[i,0,0], triangle_vertices[i,0,1]],[triangle_vertices[i,1,0], triangle_vertices[i,1,1]],[triangle_vertices[i,2,0], triangle_vertices[i,2,1]]]);
		patches.append(Polygon(v,True));
	p = PatchCollection(patches, alpha=0.4)

	p.set_array(np.array(Rdet))
	ax=plt.gca();
	ax.add_collection(p)
	plt.colorbar(p, ax=ax)
	plt.xlim([-125, -118])
	plt.ylim([36.5, 42.5])
	plt.title("Determinant of Rotation Rate Tensor")
	plt.savefig("Rdet.jpg");
	plt.close();
	
	return;


if __name__=="__main__":
	[input_file] = configure();
	[lon, lat, vel_east, vel_north, s_east, s_north] = inputs(input_file);
	[z, tri, centroids, Tdet, Rdet,Warray,Varray] = process(lon, lat, vel_east, vel_north, s_east, s_north);
	outputs(lon, lat, vel_east, vel_north, s_east, s_north, z, tri, centroids, Tdet, Rdet,Warray,Varray);


