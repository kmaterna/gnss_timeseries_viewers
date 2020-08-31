#!/usr/bin/env python
"""
Example driver for individual time series reading and writing.  
Read, Process a little, and Write Back Out.  
"""

import gps_input_pipeline
import offsets
import gps_io_functions

station="P325"
data_config_file="/Users/kmaterna/Documents/B_Research/Mendocino_Geodesy/GPS_POS_DATA/config.txt"
outfile=station+"_noearthquake.pos"

[myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(station, 'unr', data_config_file, refframe='ITRF')
newobj = offsets.remove_offsets(myData, offset_obj);   # remove antenna changes and instrument changes
newobj = offsets.remove_offsets(newobj, eq_obj);       # remove earthquakes
gps_io_functions.write_pbo_pos_file(newobj, outfile, comment="Writing a station's .pos file");
