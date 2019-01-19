#!/bin/bash
# Run the accel_mapping.gmt script
# On different areas. 

# outdir="pbo_lssq_NA/"
# name="2010"
outdir=$1
name=$2

data_file=$outdir$name.txt
./accel_map_gps.gmt $data_file -121.8 -115.0 32.2 37.6 $outdir'SoCal'$name
./accel_map_gps.gmt $data_file -125.6 -110.0 32.5 48.5 $outdir'WUS'$name
./accel_map_gps.gmt $data_file -123.5 -119.0 35.6 40.0 $outdir'SF'$name
./accel_map_gps.gmt $data_file -126.8 -121.0 38.6 43.0 $outdir'MTJ'$name
./accel_map_gps.gmt $data_file -124.6 -118.4 35.5 42.2 $outdir'NorCal'$name
./accel_map_gps.gmt $data_file -124.6 -120.4 41.2 46.2 $outdir'Oregon'$name
./accel_map_vert.gmt $data_file -125.6 -110.0 32.5 48.5 $outdir'WUS'$name
