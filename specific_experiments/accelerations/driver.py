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
Center_Coords=[-125.134, 43.829]; # Oregon


# LEAST SQUARES
# Mend_map_accel.driver(Center_Coords, "MTJ_2010", [EQ0,EQ1],[EQ1,EQ2], 'lssq'); # done
# Mend_map_accel.driver(Center_Coords, "MTJ_2010_2016", [EQ1,EQ2],[EQ3,EQ4],'lssq');  # done
# Mend_map_accel.driver(Center_Coords, "MTJ_2016", [EQ2,EQ3],[EQ3,EQ4], 'lssq'); # done
# Mend_map_accel.driver(Center_Coords, "MTJ_2014", [EQ1,EQ2],[EQ2,EQ3], 'lssq'); # done

# NOTCH FILTER
Mend_map_accel.driver(Center_Coords, "MTJ_2010", [EQ0,EQ1],[EQ1,EQ2], 'notch'); # done
Mend_map_accel.driver(Center_Coords, "MTJ_2010_2016", [EQ1,EQ2],[EQ3,EQ4],'notch');  # done
Mend_map_accel.driver(Center_Coords, "MTJ_2014", [EQ1,EQ2],[EQ2,EQ3], 'notch'); # done
Mend_map_accel.driver(Center_Coords, "MTJ_2016", [EQ2,EQ3],[EQ3,EQ4], 'notch'); # done

# GRACE-BASED
# Mend_map_accel.driver(Center_Coords, "MTJ_2010", [EQ0,EQ1],[EQ1,EQ2], 'grace'); # done
# Mend_map_accel.driver(Center_Coords, "MTJ_2014", [EQ1,EQ2],[EQ2,EQ3], 'grace'); # done

# STL FILTER
# Mend_map_accel.driver(Center_Coords, "MTJ_2010", [EQ0,EQ1],[EQ1,EQ2], 'stl'); 
# Mend_map_accel.driver(Center_Coords, "MTJ_2010_2016", [EQ1,EQ2],[EQ3,EQ4],'stl');  
# Mend_map_accel.driver(Center_Coords, "MTJ_2014", [EQ1,EQ2],[EQ2,EQ3], 'stl'); 
# Mend_map_accel.driver(Center_Coords, "MTJ_2016", [EQ2,EQ3],[EQ3,EQ4], 'stl'); 

# NO CORRECTION
# Mend_map_accel.driver(Center_Coords, "MTJ_2010", [EQ0,EQ1],[EQ1,EQ2], 'none'); 
# Mend_map_accel.driver(Center_Coords, "MTJ_2010_2016", [EQ1,EQ2],[EQ3,EQ4],'none'); 
# Mend_map_accel.driver(Center_Coords, "MTJ_2014", [EQ1,EQ2],[EQ2,EQ3], 'none'); 
# Mend_map_accel.driver(Center_Coords, "MTJ_2016", [EQ2,EQ3],[EQ3,EQ4], 'none');


# EQcoords=[-125.134, 40.829]; # The March 10, 2014 M6.8 earthquake
# EQcoords=[-122.834, 37.829]; # San Francisco Bay Area
# EQcoords=[-123.834, 39.029]; # North Bay Area
# EQcoords=[-125.134, 43.829]; # Oregon
# EQcoords=[-125.134, 46.829]; # Washington
