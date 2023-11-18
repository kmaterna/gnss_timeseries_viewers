#!/usr/bin/env python
"""
Example driver for the 'find within radius' functions
"""

from gnss_timeseries_viewers.gps_tools import load_gnss
import os

# Point to the GNSS data config file with local paths
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))  # 5 dirs up
config_file = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt');

center = [-122.0, 40.0]  # lon, lat
radius = 100  # km

database = load_gnss.create_station_repo(config_file, proc_center='pbo', refframe='NA')
stations, _ = database.search_stations_by_circle(center, radius)
print([x.name for x in stations])
