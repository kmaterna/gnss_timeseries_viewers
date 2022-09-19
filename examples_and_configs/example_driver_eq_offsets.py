# See earthquake offsets at stations
# Step 1: Determine the stations in the radius.  Then plot their time series.

import subprocess
import datetime as dt

import gps_tools.file_io.io_other
from GNSS_TimeSeries_Viewers import gps_tools as gps_tools

# CONFIG PARAMETERS FOR THIS EXPERIMENT #
params = {"data_config_file": "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt",
          "center": [-124.0, 40.15],
          "expname": 'MTJ_2014',
          "radius": 300,  # in km
          "proc_center": 'cwu',
          "refframe": 'NA',
          "blacklist": ()};


def driver():
    myparams, stations = configure();
    [data, _, eq_list, _] = gps_tools.gps_input_pipeline.multi_station_inputs(stations, myparams["blacklist"],
                                                                              myparams["proc_center"],
                                                                              myparams["refframe"],
                                                                              myparams["data_config_file"],
                                                                              myparams["distances"]);
    # FOR 2014 earthquake
    station_vectors = gps_tools.offsets.offset_to_vel_object(eq_list, data, myparams["refframe"], myparams["proc_center"],
                                                             target_date=dt.datetime.strptime("20140310", "%Y%m%d"));
    gps_tools.gps_vel_pygmt_plots.simple_pygmt_plot(station_vectors, myparams["outdir"] + "/mtj_vector_map.png",
                                                    symsize=0.1, vector_scale_info=(0.5, 10, "10 mm"));
    gps_tools.file_io.io_other.write_stationvel_file(station_vectors, myparams["outdir"] + '/mtj_vectors.txt',
                                                     metadata='cwu_NA look-up table for 03-10-2014');
    return;


def configure():
    stations, _, _, distances = gps_tools.stations_within_radius.get_stations_within_radius(params["data_config_file"],
                                                                                            params["center"],
                                                                                            params["radius"],
                                                                                            network=params["proc_center"]);
    outdir = params["expname"] + "_" + params["proc_center"]
    subprocess.call(["mkdir", "-p", outdir], shell=False);
    params["outdir"] = outdir;
    outname = params["expname"] + "_" + str(params["center"][0]) + "_" + str(params["center"][1]) + "_" + str(
        params["radius"])
    params["outname"] = outname;
    return params, stations;


if __name__ == "__main__":
    driver();
