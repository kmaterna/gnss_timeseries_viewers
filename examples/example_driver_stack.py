#!/usr/bin/env python
"""
Example driver for stacked time series viewing. 
"""

import os
from gnss_timeseries_viewers.gps_tools import gps_stack

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # 4 dirs up
data_config_file = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt')
expname = 'NBay'
center = [-122.0, 38.0]  # lon, lat
radius = 30  # km
proc_center = 'cwu'  # Which datastream do you want?
refframe = 'NA'  # Which reference frame?
outdir = expname + "_" + proc_center
gps_stack.driver(data_config_file, expname, center, radius, proc_center, refframe, outdir)
