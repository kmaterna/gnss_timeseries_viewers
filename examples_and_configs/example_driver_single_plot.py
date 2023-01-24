#!/usr/bin/env python
"""
Example driver for individual time series viewing. 
"""

import sys
from GNSS_TimeSeries_Viewers.gps_tools import single_station_tsplot

station = "P325"
data_config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"

if len(sys.argv) >= 2:
    station = sys.argv[1];  # you can type in the name of a station in the run string (if you want)

single_station_tsplot.view_single_station(station, data_config_file=data_config_file, offsets_remove=1,
                                          earthquakes_remove=0, outliers_remove=1, seasonals_remove=1, outliers_def=15,
                                          seasonals_type='lssq', datasource='cwu', refframe='NA', outdir="");
