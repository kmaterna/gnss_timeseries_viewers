"""
Download GNSS time series from the University of Nevada Reno (http://geodesy.unr.edu/PlugNPlayPortal.php).
Inputs are based on a cache file (for coordinates) and a geographic region.
Will download in several reference frames, such as NA and IGS14.
Call this script from the directory where time series will live.
Example runstring:
Script by K. Materna, 2020
"""

import numpy as np
import sys
from subprocess import call

def configure(argv):
    """
    Before 2018: ref_frame is either ISG08 or NA12
    After 2018: ref_frame is either IGS14 or NA
    For most of California, it's about 1200 stations and about 500Mb. It takes ~10 minutes.
    For the western US, it's about 2600 stations and ~1 hour.
    """
    if len(argv) != 2:
        print("\n\nError! Please provide get_unr_time_series.py with the name of the coordinate cache file. ")
        print("ex: get_unr_time_series.py ../systems/path/to/UNR_coords.txt\n")
        sys.exit(0);
    else:
        coordfile = argv[1];
    ref_frame1 = "NA";
    ref_frame2 = "IGS14";
    latlon_box = [-125.0, -110, 32.0, 49.0];  # A LARGE BOX INCLUDING ALL OF WUS
    return coordfile, ref_frame1, ref_frame2, latlon_box;

def get_stations(coordfile, latlon_box):
    station_names = [];
    data = np.genfromtxt(coordfile, skip_header=2, usecols=(0, 1, 2), dtype=None);
    for item in data:
        lat = item[1];
        lon = item[2];
        if lon > 180:
            lon = lon-360;
        if latlon_box[0] <= lon <= latlon_box[1]:
            if latlon_box[2] <= lat <= latlon_box[3]:
                onestation = item[0].decode('UTF-8');
                station_names.append(onestation);
    print("Returning %d stations for downloading. " % len(station_names));
    return station_names;

def download_stations(stations, ref_frame):
    if ref_frame == "NA":
        ref_frame_local = "NA";
    else:
        ref_frame_local = ref_frame;
    for mystation in stations:
        print("Downloading station %s " % mystation);
        if ref_frame == "IGS14":
            call('wget http://geodesy.unr.edu/gps_timeseries/tenv3/'+ref_frame+'/'+mystation+'.tenv3 -O ' +
                 mystation+'.'+ref_frame+'.tenv3', shell=True);
        else:
            call('wget http://geodesy.unr.edu/gps_timeseries/tenv3/'+'plates'+'/'+ref_frame+'/' +
                 mystation+'.'+ref_frame+'.tenv3 -O '+mystation+'.'+ref_frame_local+'.tenv3', shell=True);
    return;


if __name__ == "__main__":
    coordfile, ref_frame1, ref_frame2, latlon_box = configure(sys.argv);
    stations = get_stations(coordfile, latlon_box);
    download_stations(stations, ref_frame1);  # NA first
    download_stations(stations, ref_frame2);  # IGS14 next
