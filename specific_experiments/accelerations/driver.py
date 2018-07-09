"""
Driver for making maps of GPS accelerations before/after a certain time. 
In this case, EQcoords is just the center of the map. 
"""

import Mend_map_accel


basename="_accel_3yrs_2011.ps"

# EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
# Mend_map_accel.driver(EQcoords, "MTJ"+basename);

# EQcoords=[-122.834, 37.829]; # San Francisco Bay Area
# Mend_map_accel.driver(EQcoords, "SF"+basename);

EQcoords=[-123.834, 39.029]; # North Bay Area
Mend_map_accel.driver(EQcoords, "NB"+basename);

# EQcoords=[-125.134, 43.829]; # Oregon
# Mend_map_accel.driver(EQcoords, "OR"+basename);

# # EQcoords=[-125.134, 46.829]; # Washington
# # Mend_map_accel.driver(EQcoords, "WA"+basename);

# EQcoords=[-120.134, 33.829]; # SoCal
# Mend_map_accel.driver(EQcoords, "SoCal"+basename);
