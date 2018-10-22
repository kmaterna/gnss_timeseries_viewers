"""
Driver for making maps of GPS accelerations before/after a certain time. 
In this case, EQcoords is just the center of the map. 
"""

import Mend_map_accel

EQ0="20050101"
EQ1="20100110"
EQ2="20140310"
EQ3="20161208"
EQ4="20180905"

# Center_Coords=[-125.134, 40.829]; # MTJ
# Center_Coords=[-122.834, 37.829]; # San Francisco Bay Area
# Center_Coords=[-125.134, 43.829]; # Oregon
Center_Coords=[-124.0, 38.5]; # Western US everything
filter_type='lssq'; # OPTIONS: lssq, notch, stl, grace, none. 
size='huge';  # OPTIONS: small, medium, huge
network='unr';


# LEAST SQUARES
# Mend_map_accel.driver(Center_Coords, "MTJ_2010", [EQ0,EQ1],[EQ1,EQ2], filter_type, size, network); # done
# Mend_map_accel.driver(Center_Coords, "MTJ_2010_2016", [EQ1,EQ2],[EQ3,EQ4],filter_type, size, network);  # done
# Mend_map_accel.driver(Center_Coords, "MTJ_2014", [EQ1,EQ2],[EQ2,EQ3], filter_type, size, network); # done
Mend_map_accel.driver(Center_Coords, "MTJ_2016", [EQ2,EQ3],[EQ3,EQ4], filter_type, size, network); # done




