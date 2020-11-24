"""
Driver for individual time series or stacked time series viewing. 
"""

import sys
import single_station_tsplot

station="P349"

if len(sys.argv) >=2:
	station=sys.argv[1];  # you can type in the name of a station instead (if you want)

single_station_tsplot.view_single_station(station, 
	offsets_remove=1, earthquakes_remove=1, 
	outliers_remove=1, seasonals_remove=1, outliers_def=15, seasonals_type='grace', datasource='unr', refframe='NA', 
	data_config_file="/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt");


# 9/27/2018: NOTES
# P319 has a strange north velocity. 
# Grep P323 didn't return nice things. 
# BRAW has problems with a gap. 
# starttime=dt.datetime.strptime("20090505","%Y%m%d") used for BRAW experiment
