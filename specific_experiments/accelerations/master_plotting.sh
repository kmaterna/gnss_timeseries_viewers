# Run the accel_mapping.gmt script
# On different areas. 

# outdir="NAPA_lssq_UNR/"
# data_file=$outdir"NAPA_lssq.txt"
# out_base="_NAPA_lssq"

outdir="Whole_WUS_lssq_UNR_itrf/"
data_file=$outdir"MTJ_2016_lssq.txt"
out_base="_2016_lssq"

./accel_map_gps.gmt $data_file -121.8 -115.0 32.2 37.6 $outdir'SoCal'$out_base
./accel_map_gps.gmt $data_file -125.6 -110.0 32.5 48.5 $outdir'WUS'$out_base
./accel_map_gps.gmt $data_file -123.5 -119.0 35.6 40.0 $outdir'SF'$out_base
./accel_map_gps.gmt $data_file -125.2 -121.0 38.6 43.0 $outdir'MTJ'$out_base
./accel_map_gps.gmt $data_file -124.6 -118.4 35.5 42.2 $outdir'NorCal'$out_base
./accel_map_gps.gmt $data_file -124.6 -120.4 41.2 46.2 $outdir'Oregon'$out_base