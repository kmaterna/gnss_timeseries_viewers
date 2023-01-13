#!/usr/bin/env python
"""
Example driver for individual time series reading and writing.  
Read, Process a little, and Write Back Out.  
"""
import gps_tools.file_io.io_nota
from GNSS_TimeSeries_Viewers.gps_tools import offsets, load_gnss

station_name = "P325"
data_config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"
outfile = station_name + "_noearthquake.pos"

database = load_gnss.create_station_repo(data_config_file, 'ITRF', 'unr')
[myData, offset_obj, eq_obj] = database.load_station(station_name);
newobj = offsets.remove_offsets(myData, offset_obj);   # remove antenna changes and instrument changes
newobj = offsets.remove_offsets(newobj, eq_obj);       # remove earthquakes
gps_tools.file_io.io_nota.write_pbo_pos_file(newobj, outfile, comment="Writing a station's .pos file");
