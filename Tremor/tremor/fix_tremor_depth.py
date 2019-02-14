# Goal: 
# Read McCrory et al. (2012) geometry
# Use this surface to estimate the depths of all the tremor
# Based on their lat/lon

import numpy as np
import matplotlib.pyplot as plt 
import subprocess
import netcdf_read_write
import tremor_io
import tremor_tools


def read_csz_model():
	grdname="../slab_geometry/cas_slab1.0_clip.grd";
	[xdata, ydata, zdata] = netcdf_read_write.read_grd_xyz(grdname);
	xdata=np.flipud(xdata);
	ydata=np.flipud(ydata);
	zdata=np.flipud(zdata);
	return xdata, ydata, zdata;

def compute_depths(tremor, xdata, ydata, zdata):
	coords=[];
	for i in range(len(tremor.lonarray)):
		coords.append([tremor.lonarray[i], tremor.latarray[i]]);
	depths=tremor_tools.get_depth_projection(coords,xdata,ydata,zdata);	
	tremor_with_depths=tremor_tools.associate_depths(tremor, depths);
	return tremor_with_depths;

def make_plots(xdata, ydata, zdata, tremor):
	plt.contourf(xdata,ydata,zdata,30);
	plt.colorbar();
	plt.scatter(tremor.lonarray,tremor.latarray,0.5,c=tremor.depth);
	plt.savefig('depth_test.eps');



if __name__=="__main__":
	tremor = tremor_io.read_input_tremor("ide");
	xdata, ydata, zdata = read_csz_model();
	tremor_with_depths = compute_depths(tremor, xdata, ydata, zdata);
	make_plots(xdata, ydata, zdata, tremor_with_depths);



