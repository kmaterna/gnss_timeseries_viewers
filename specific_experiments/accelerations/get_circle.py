# Python script to get the circle of points at a given radius surrounding a coordinate

import numpy as np 
import add_coords_km


[lon0, lat0] = [-124.81, 40.53];
theta=np.arange(0,2*np.pi,0.01);
radius=350; 
ofile=open(str(radius)+"rad.txt",'w');

for i in range(len(theta)):
	return_coords = add_coords_km.add_vector_to_coords(lon0, lat0, radius*np.cos(theta[i]), radius*np.sin(theta[i]));
	ofile.write('%f %f\n' % (return_coords[0], return_coords[1]) );
ofile.close();
