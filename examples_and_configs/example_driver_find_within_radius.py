#!/usr/bin/env python
"""
Example driver for the 'find within radius' functions
"""

from GNSS_TimeSeries_Viewers.gps_tools import load_gnss

config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"
center = [-122.0, 40.0]  # lon, lat
radius = 100;  # km

database = load_gnss.create_station_repo(config_file, proc_center='pbo', refframe='NA');
stations, _ = database.search_stations_by_circle(center, radius);
print([x.name for x in stations]);
