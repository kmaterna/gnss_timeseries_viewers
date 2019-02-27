# Goal: 
# Read McCrory et al. (2012) geometry
# Use this surface to estimate the depths of all the tremor based on their lat/lon
# Run the more complicated plots using the depths instead of simple boxes. 

import numpy as np
import matplotlib.pyplot as plt 
import datetime as dt 
import subprocess
import sys
import netcdf_read_write
import tremor_io
import tremor_tools
import tremor_plots

def read_csz_model():
	grdname="../slab_geometry/cas_slab1.0_clip.grd";
	[xdata, ydata, zdata] = netcdf_read_write.read_grd_xyz(grdname);
	xdata=np.flipud(xdata);
	ydata=np.flipud(ydata);
	zdata=np.flipud(zdata);
	return xdata, ydata, zdata;

def combine_tremor(tremor_type):
	if tremor_type=="wech_custom":
		tremor_website = tremor_io.read_input_tremor("wech");
		tremor_later = tremor_io.read_input_tremor("wech_custom");
		transition_time_1 = dt.datetime.strptime("2015-01-01","%Y-%m-%d");
		transition_time_2 = dt.datetime.strptime("2018-01-01","%Y-%m-%d");
		box_interest=[-125,-120,38,43];
		tremor_website1 = tremor_tools.restrict_to_box(tremor_website, box_interest, tremor_website.dtarray[0], transition_time_1);
		tremor_website2 = tremor_tools.restrict_to_box(tremor_website, box_interest, transition_time_2, tremor_website.dtarray[-1]);
		tremor_total = tremor_tools.concatonate_tremor(tremor_website1, tremor_later);
		tremor_total = tremor_tools.concatonate_tremor(tremor_total, tremor_website2)
	else:
		tremor_total = tremor_io.read_input_tremor(tremor_type);
	return tremor_total;

def compute_depths(tremor, xdata, ydata, zdata):
	coords=[];
	for i in range(len(tremor.lonarray)):
		coords.append([tremor.lonarray[i], tremor.latarray[i]]);
	depths=tremor_tools.get_depth_projection(coords,xdata,ydata,zdata);	
	tremor_with_depths=tremor_tools.associate_depths(tremor, depths);
	return tremor_with_depths;

def make_plots(xdata, ydata, zdata, tremor):
	# Make plot of the csz model vs. tremor depths (IT WORKS!)
	plt.contourf(xdata,ydata,zdata,30);
	plt.colorbar();
	plt.scatter(tremor.lonarray,tremor.latarray,0.5,c=tremor.depth);
	plt.savefig('depth_test.eps');


if __name__=="__main__":
	tremor_type="wech_custom";
	tremor = combine_tremor(tremor_type);
	xdata, ydata, zdata = read_csz_model();
	tremor_with_depths = compute_depths(tremor, xdata, ydata, zdata);
	tremor_plots.complex_plot_depths(tremor_with_depths,tremor_type);
	# After this, you must go and make the GMT plots of the tremor (tremor_depth_ranges.sh)
