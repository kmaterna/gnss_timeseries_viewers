#!/usr/bin/env python
"""
Driver for individual time series reading and writing. 
"""

import sys

station="P325"

if len(sys.argv) >=2:
	station=sys.argv[1];  # you can type in the name of a station instead (if you want)

