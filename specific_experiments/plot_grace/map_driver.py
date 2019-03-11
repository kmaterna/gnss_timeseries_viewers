# Take the GRACE time series and make velocity change maps. 

import subprocess
import Mend_map_accel

grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"
out_dir="maps/";
EQ0="20050101"
EQ1="20100110"
EQ2="20140310"
EQ3="20161208"
EQ4="20180905"

EQ1="20100110"
EQ2="20140310"
EQ3="20161210"

dataset="MTJ_2014"
Mend_map_accel.grace_driver([EQ1,EQ2],[EQ2,EQ3], grace_dir, dataset, out_dir);
subprocess.call('../accelerations/accel_map_grace.gmt '+out_dir+dataset+'.txt '+' -126.8 -119.4 35.9 44.2 '+out_dir+'/NorCal_'+dataset,shell=True);

# dataset="MTJ_2010"
# Mend_map_accel.grace_driver([EQ0,EQ1],[EQ1,EQ2], grace_dir, dataset, out_dir);
# subprocess.call('../accelerations/accel_map_grace.gmt '+out_dir+dataset+'.txt '+' -126.8 -119.4 37.8 42.2 '+out_dir+'/NorCal_'+dataset,shell=True);

