#!/usr/bin/env python
"""
Example driver for stacked time series viewing. 
"""

import datetime as dt
from GNSS_TimeSeries_Viewers.gps_tools import gps_stack

data_config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"
expname = 'NBay';
center = [-122.0, 38.0];  # lon, lat
radius = 40;  # km
proc_center = 'cwu';  # Which datastream do you want? 
refframe = 'NA';  # Which reference frame? 
outdir = expname + "_" + proc_center
must_include_window = [dt.datetime.strptime("20100310", "%Y%m%d"), dt.datetime.strptime("20140310", "%Y%m%d")];
# time window we must include
gps_stack.driver(data_config_file, expname, center, radius, proc_center, refframe, outdir, must_include_window);
