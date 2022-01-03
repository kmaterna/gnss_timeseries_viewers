#!/usr/bin/env python
"""
Example driver for the 'find within radius' functions
"""

from GNSS_TimeSeries_Viewers.gps_tools import stations_within_radius

data_config_file="/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"
center = [-122.0, 40.0]  # lon, lat
radius = 100;  # km

stations, lon, lat, distances = stations_within_radius.get_stations_within_radius(data_config_file, center, radius, network='pbo');
print(stations);
