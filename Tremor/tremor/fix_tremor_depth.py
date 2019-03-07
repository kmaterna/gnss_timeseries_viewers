# Goal: 
# Read McCrory et al. (2012) geometry
# Use this surface to estimate the depths of all the tremor based on their lat/lon
# Run the more complicated plots using the depths instead of simple boxes. 

import numpy as np
import matplotlib.pyplot as plt 
import datetime as dt 
import subprocess
import sys
import tremor_tools
import tremor_plots




if __name__=="__main__":
	tremor_type="wech_custom";
	tremor_with_depths = tremor_tools.read_custom_tremor(tremor_type);
	tremor_plots.complex_plot_depths(tremor_with_depths,tremor_type);
	sys.exit(0);
	# After this, you must go and make the GMT plots of the tremor (tremor_depth_ranges.sh)

	depth_interest=[10, 65];
	box_interest = [-125, -121, 40.2, 40.8];
	# dt1_start = dt.datetime.strptime("20140222","%Y%m%d");
	dt1_start = dt.datetime.strptime("20140101","%Y%m%d");
	dt1_end = dt.datetime.strptime("20140303","%Y%m%d");
	dt2_start = dt.datetime.strptime("20140905","%Y%m%d");
	dt2_end = dt.datetime.strptime("20140918","%Y%m%d");
	dt3_start = dt.datetime.strptime("20150423","%Y%m%d");
	dt3_end = dt.datetime.strptime("20150517","%Y%m%d");	
	dt4_start = dt.datetime.strptime("20151221","%Y%m%d");
	dt4_end = dt.datetime.strptime("20151230","%Y%m%d");
	dt5_start = dt.datetime.strptime("20160810","%Y%m%d");
	dt5_end = dt.datetime.strptime("20160820","%Y%m%d");
	dt6_start = dt.datetime.strptime("20170422","%Y%m%d");
	dt6_end = dt.datetime.strptime("20170502","%Y%m%d");
	dt7_start = dt.datetime.strptime("20171221","%Y%m%d");
	dt7_end = dt.datetime.strptime("20180105","%Y%m%d");
	interval_list=[[dt1_start, dt1_end], [dt2_start, dt2_end], [dt3_start, dt3_end], [dt4_start, dt4_end], [dt5_start, dt5_end], [dt6_start, dt6_end], [dt7_start, dt7_end]];
	tremor_plots.histogram_depths(tremor_with_depths, interval_list, box_interest, depth_interest);

	# Get intervals. 
	# Compute histogram for each interval. 

