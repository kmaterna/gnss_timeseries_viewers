#!/usr/bin/env python
"""
Example driver for stacked time series viewing. 
"""

from GNSS_TimeSeries_Viewers.gps_tools import gps_stack

data_config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"
expname = 'NBay';
center = [-122.0, 38.0];  # lon, lat
radius = 30;  # km
proc_center = 'cwu';  # Which datastream do you want? 
refframe = 'NA';  # Which reference frame? 
outdir = expname + "_" + proc_center
gps_stack.driver(data_config_file, expname, center, radius, proc_center, refframe, outdir);
