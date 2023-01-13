#!/usr/bin/python3

# See earthquake offsets at stations
# Determine the stations in the radius.  Then identify their offsets.
# This could be further automated with a strongly typed params object

import subprocess
import datetime as dt
from GNSS_TimeSeries_Viewers import gps_tools as gt
import GNSS_TimeSeries_Viewers.gps_tools.file_io.io_other as io_other

# CONFIG PARAMETERS FOR THIS EXPERIMENT #
params = {"data_config_file": "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt",
          "center": [-124.0, 41.15],
          "expname": 'MTJ_2014',
          "radius": 200,  # in km
          "proc_center": 'cwu',
          "refframe": 'NA',
          "blacklist": ()};


def driver():
    myparams, database, stations = configure();
    [data, _, eq_list] = database.load_stations(stations);
    # FOR 2014 earthquake
    station_vectors = gt.offsets.offset_to_vel_object(eq_list, data, myparams["refframe"], myparams["proc_center"],
                                                      target_date=dt.datetime.strptime("20140310", "%Y%m%d"));
    gt.pygmt_plots.simple_pygmt_plot(station_vectors, myparams["outdir"] + "/mtj_vector_map.png",
                                     symsize=0.1, vector_scale_info=(0.5, 10, "10 mm"));
    io_other.write_stationvel_file(station_vectors, myparams["outdir"] + '/mtj_vectors.txt',
                                   metadata='cwu_NA look-up table for 03-10-2014');
    return;


def configure():
    database = gt.load_gnss.create_station_repo(params['data_config_file'], proc_center=params["proc_center"],
                                                refframe=params['refframe']);
    stations, _ = database.search_stations_by_circle(params['center'], params['radius']);
    outdir = params["expname"] + "_" + params["proc_center"]
    subprocess.call(["mkdir", "-p", outdir], shell=False);
    params["outdir"] = outdir;
    return params, database, [x.name for x in stations];


if __name__ == "__main__":
    driver();
