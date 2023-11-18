#!/usr/bin/env python

# See earthquake offsets at stations
# Determine the stations in the radius.  Then identify their offsets.

import os
import datetime as dt
from gnss_timeseries_viewers.gps_tools import offsets, pygmt_plots, load_gnss
import gnss_timeseries_viewers.gps_tools.file_io.io_other as io_other

# CONFIG PARAMETERS FOR THIS EXPERIMENT #
# Point to the GNSS data config file with local paths
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # 4 dirs up
data_config_file = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt')
params = {"center": [-124.0, 41.15],
          "expname": 'MTJ_2014',
          "radius": 200,  # in km
          "proc_center": 'cwu',
          "refframe": 'NA',
          "blacklist": ()}


def driver():
    myparams, database, stations = configure()
    [data, _, eq_list] = database.load_stations(stations)
    offsetpts = offsets.table_offset_to_velfield(data, eq_list, target_date=dt.datetime.strptime("20140310", "%Y%m%d"))
    pygmt_plots.simple_pygmt_plot(offsetpts, myparams["outdir"] + "/mtj_vector_map.png",
                                  symsize=0.1, vector_scale_info=(0.5, 10, "10 mm"))
    io_other.write_stationvel_file(offsetpts, myparams["outdir"] + '/mtj_vectors.txt',
                                   metadata='cwu_NA look-up table for 03-10-2014')
    return


def configure():
    database = load_gnss.create_station_repo(data_config_file, proc_center=params["proc_center"],
                                             refframe=params['refframe'])
    stations, _ = database.search_stations_by_circle(params['center'], params['radius'])
    outdir = params["expname"] + "_" + params["proc_center"]
    os.makedirs(outdir, exist_ok=True)
    params["outdir"] = outdir
    return params, database, [x.name for x in stations]


if __name__ == "__main__":
    driver()
