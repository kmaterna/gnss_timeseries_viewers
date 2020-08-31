#!/usr/bin/env python
"""
Example driver for individual time series viewing. 
"""

import sys
import single_station_tsplot

station="P325"

if len(sys.argv) >=2:
	station=sys.argv[1];  # you can type in the name of a station in the run string (if you want)

single_station_tsplot.view_single_station(station, 
	offsets_remove=1, earthquakes_remove=0, 
	outliers_remove=1, seasonals_remove=1, outliers_def=15, seasonals_type='lssq', datasource='cwu', refframe='NA', 
	data_config_file="/Users/kmaterna/Documents/B_Research/Mendocino_Geodesy/GPS_POS_DATA/config.txt", outdir="");

