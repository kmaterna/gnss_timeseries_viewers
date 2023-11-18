#!/usr/bin/env python
"""
Example driver for individual time series reading and writing.  
Read, Process a little, and Write Back Out.  
"""
from gnss_timeseries_viewers.gps_tools.file_io import io_nota
from gnss_timeseries_viewers.gps_tools import offsets, load_gnss
import os

station_name = "P325"
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # 4 dirs up
data_config_file = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt')
outfile = station_name + "_noearthquake.pos"

database = load_gnss.create_station_repo(data_config_file, 'ITRF', 'unr')
[myData, offset_obj, eq_obj] = database.load_station(station_name)
newobj = offsets.remove_offsets(myData, offset_obj)   # remove antenna changes and instrument changes
newobj = offsets.remove_offsets(newobj, eq_obj)       # remove earthquakes
io_nota.write_pbo_pos_file(newobj, outfile, comment="Writing a station's .pos file")
