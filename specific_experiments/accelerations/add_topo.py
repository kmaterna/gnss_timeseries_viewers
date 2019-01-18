# Add topography

import subprocess

outdir="Whole_WUS_lssq_PBO_nam/";
dataset="2014_lssq";



datafile=outdir+"MTJ_"+dataset+".txt";
subprocess.call('./accel_map_gps.gmt '+datafile+' -126.7 -121.0 38.8 42.9 '+outdir+'/MTJ_'+dataset,shell=True);  
 # MTJ: -125.5 -121.0 38.8 42.9
subprocess.call('./accel_map_gps.gmt '+datafile+' -124.6 -118.4 35.5 42.2 '+outdir+'/NorCal_'+dataset,shell=True);
subprocess.call('./accel_map_gps.gmt '+datafile+' -121.8 -115.0 32.2 37.6 '+outdir+'/SoCal_'+dataset,shell=True);
subprocess.call('./accel_map_gps.gmt '+datafile+' -125.6 -110.0 32.5 48.5 '+outdir+'/WUS_'+dataset,shell=True);
subprocess.call('./accel_map_gps.gmt '+datafile+' -123.5 -119.0 35.6 40.0 '+outdir+'/SF_'+dataset,shell=True);
subprocess.call('./accel_map_gps.gmt '+datafile+' -124.8 -120.4 41.2 46.2 '+outdir+'/Oregon_'+dataset,shell=True);