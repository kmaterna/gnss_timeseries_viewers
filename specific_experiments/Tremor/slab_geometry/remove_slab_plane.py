# Goal: 
# Read McCrory et al. (2012) geometry
# Estimate a best-fitting plane
# Subtract plane from geometry
# Plot the residuals.

import numpy as np
import matplotlib.pyplot as plt 
import subprocess
import netcdf_read_write

grdname="cas_slab1.0_clip.grd";
xyzname="casc_geometry.xyz";
detrend_xyz="detrended_geometry.xyz";
detrend_grd="detrended_N4.grd";

top_range=47;  # degrees north. 

[xdata, ydata, zdata] = netcdf_read_write.read_grd_xyz(grdname);
xdata=np.flipud(xdata);

# plt.figure();
# plt.contourf(xdata, ydata, zdata);
# plt.colorbar();
# plt.savefig('test.eps');

outfile=open(xyzname,'w');
for i in range(len(ydata)):
	for j in range(len(xdata)):
		if ydata[i]<top_range:
			if ~ np.isnan(zdata[i][j]):
				outfile.write("%f %f %f\n" % (xdata[j], ydata[i], zdata[i][j]) );
outfile.close();

# Remove a plane
order=str(4);
subprocess.call(["gmt","trend2d",xyzname,"-Fxyr","-N"+order,"-V",">",detrend_xyz]);
subprocess.call(["gmt xyz2grd "+detrend_xyz+" -G"+detrend_grd+" `gmt grdinfo -I "+grdname+"` -R-128.5/-120.5/39/"+str(top_range)],shell=True);
subprocess.call(["./topo_map1.gmt"],shell=False);

# [y,x,r]=np.loadtxt(detrend_name, unpack=True);
# plt.figure();
# plt.scatter(x,y,c=r,s=1);
# plt.colorbar();
# plt.savefig('resids.eps');