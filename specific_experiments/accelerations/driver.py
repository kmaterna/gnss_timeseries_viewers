"""
Driver for making maps of GPS accelerations before/after a certain time. 
In this case, EQcoords is just the center of the map. 
"""

import Mend_map_accel


EQ0="20050101"
EQ1="20100110"
EQ2="20140314"
EQ3="20161208"
EQ4="20180505"

# EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
# Mend_map_accel.driver(EQcoords, "MTJ"+basename);

# EQcoords=[-122.834, 37.829]; # San Francisco Bay Area
# Mend_map_accel.driver(EQcoords, "SF"+basename);

# EQcoords=[-123.834, 39.029]; # North Bay Area
# Mend_map_accel.driver(EQcoords, "NB"+basename);

# EQcoords=[-125.134, 43.829]; # Oregon
# Mend_map_accel.driver(EQcoords, "OR"+basename);

# # EQcoords=[-125.134, 46.829]; # Washington
# # Mend_map_accel.driver(EQcoords, "WA"+basename);

EQcoords=[-125.134, 40.829]; # MTJ
Mend_map_accel.driver(EQcoords, "MTJ_2010.ps", [EQ0,EQ1],[EQ1,EQ2]);
Mend_map_accel.driver(EQcoords, "MTJ_2014.ps", [EQ1,EQ2],[EQ2,EQ3]);
Mend_map_accel.driver(EQcoords, "MTJ_2016.ps", [EQ2,EQ3],[EQ3,EQ4]);
