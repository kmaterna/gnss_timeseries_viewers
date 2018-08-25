"""
Driver for making maps of GPS accelerations before/after a certain time. 
In this case, EQcoords is just the center of the map. 
"""

import Mend_map_accel


EQ0="20050101"
EQ1="20100110"
EQ2="20140314"
EQ3="20161208"
EQ4="20180905"


Center_Coords=[-125.134, 40.829]; # MTJ
# Mend_map_accel.driver(Center_Coords, "MTJ_2010.ps", [EQ0,EQ1],[EQ1,EQ2],'horizontal');
# Mend_map_accel.driver(Center_Coords, "MTJ_2014.ps", [EQ1,EQ2],[EQ2,EQ3],'horizontal');
# Mend_map_accel.driver(Center_Coords, "MTJ_2016.ps", [EQ2,EQ3],[EQ3,EQ4],'horizontal');
# Mend_map_accel.driver(Center_Coords, "MTJ_2010_2016.ps", [EQ1,EQ2],[EQ3,EQ4],'horizontal');

Mend_map_accel.driver(Center_Coords, "MTJ_2016_vert.ps", [EQ2,EQ3],[EQ3,EQ4],'vertical');
Mend_map_accel.driver(Center_Coords, "MTJ_2014_vert.ps", [EQ1,EQ2],[EQ2,EQ3],'vertical');


# Notes: 
# Trying by multiplying the verticals by 0.3

# EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
# EQcoords=[-122.834, 37.829]; # San Francisco Bay Area
# EQcoords=[-123.834, 39.029]; # North Bay Area
# EQcoords=[-125.134, 43.829]; # Oregon
# # EQcoords=[-125.134, 46.829]; # Washington
