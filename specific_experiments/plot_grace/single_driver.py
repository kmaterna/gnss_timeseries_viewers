


import grace_ts_functions

station_name='MUSB';
grace_dir="../../GPS_POS_DATA/GRACE_loading_model/"
filename=grace_dir+"scaled_"+station_name+"_PREM_model_ts.txt";
out_dir="stations/";


grace_ts_functions.plot_grace(station_name, filename, out_dir);
