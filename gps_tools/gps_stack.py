"""
Viewing a stack of stations
  Step 1: Determine the stations in the radius.
  Step 2: Read them in. Make a list of timeseries objects in one large dataobject.
  Step 3: Compute: Remove outliers, earthquakes, and eventually trend from the data.
  Step 4: Plot in order of increasing latitude, colored by how close they are to the central point
"""

import subprocess
from . import gps_seasonal_removals, offsets, load_gnss, vel_functions, pygmt_plots
from . import outputs_gps_stacks as out_stack
from .file_io import config_io, io_other
import datetime as dt


def driver(data_config_file, expname, center, radius, proc_center, refframe, outdir):
    outname = configure(expname, center, radius, outdir);

    database, stations, distances = build_database(data_config_file, proc_center, refframe, center, radius);
    [dataobj_list, offsetobj_list, eqobj_list] = database.load_stations(stations);
    [_, _, no_offsets_no_trends, no_offsets_no_trends_no_seasons,
     sorted_distances] = compute(dataobj_list, offsetobj_list, eqobj_list, distances, data_config_file);
    start = dt.datetime.strptime("20050505", "%Y%m%d")

    # A series of output options that can be chained together, selected or unselected, etc.
    out_stack.horizontal_full_ts(no_offsets_no_trends, sorted_distances, outname=outdir+"/"+outname+'_TS_noeq.png',
                                 start_time_plot=start);
    out_stack.horizontal_full_ts(no_offsets_no_trends_no_seasons, sorted_distances,
                                 outname=outdir + "/" + outname + '_TS_noeq_noseasons.png', start_time_plot=start);
    out_stack.vertical_full_ts(no_offsets_no_trends_no_seasons, sorted_distances,
                               outname=outdir + "/" + outname + '_TS_vertical.png', start_time_plot=start);

    out_stack.horizontal_filtered_plots(no_offsets_no_trends_no_seasons, sorted_distances,
                                        outname=outdir + "/" + outname + '_TS_horiz_filt.png', start_time_plot=start);
    out_stack.vertical_filtered_plots(no_offsets_no_trends_no_seasons, sorted_distances,
                                      outname=outdir + "/" + outname + '_TS_vert_detrended_filt.png',
                                      start_time_plot=start);

    pygmt_plots.map_ts_objects(dataobj_list, outdir+"/"+outname+'_map.png', center=center);
    return;


def build_database(data_config_file, proc_center, refframe, center, radius):
    # Set up the stacking process
    database = load_gnss.create_station_repo(data_config_file, proc_center=proc_center, refframe=refframe);
    stations, distances = database.search_stations_by_circle(center, radius, basic_clean=True);
    data_config = config_io.read_config_file(data_config_file);
    blacklist = io_other.read_blacklist(data_config["blacklist"]);
    stations, distances = vel_functions.remove_blacklist_paired_data(stations, blacklist, matching_data=distances);
    return database, [x.name for x in stations], distances;


def configure(expname, center, radius, outdir):
    subprocess.call(["mkdir", "-p", outdir], shell=False);
    outname = expname + "_" + str(center[0]) + "_" + str(center[1]) + "_" + str(radius)
    return outname;


def compute(dataobj_list, offsetobj_list, eqobj_list, distances, data_config_file):
    latitudes_list = [i.coords[1] for i in dataobj_list];
    sorted_objects = [x for _, x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data.
    sorted_offsets = [x for _, x in sorted(zip(latitudes_list, offsetobj_list))];  # the raw, sorted data.
    sorted_eqs = [x for _, x in sorted(zip(latitudes_list, eqobj_list))];  # the raw, sorted data.
    sorted_distances = [x for _, x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

    detrended_objects, no_offset_objects = [], [];
    no_offsets_no_trends, no_offsets_no_trends_no_seasons = [], [];

    # Detrended objects (or objects with trends and no offsets; depends on what you want.)
    for i in range(len(sorted_objects)):
        newobj = gps_seasonal_removals.make_detrended_ts(sorted_objects[i], 0, 'lssq', data_config_file);
        detrended_objects.append(newobj);  # still has offsets, doesn't have trends

        newobj = offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
        newobj = offsets.remove_offsets(newobj, sorted_eqs[i]);
        no_offset_objects.append(newobj);  # still has trends, doesn't have offsets

    # # Objects with no earthquakes or seasonals
    for i in range(len(sorted_objects)):
        # Remove the steps earthquakes
        newobj = offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
        newobj = offsets.remove_offsets(newobj, sorted_eqs[i]);
        newobj = newobj.remove_outliers(20);  # 20mm outlier definition

        # The detrended TS without earthquakes
        stage1obj = gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq', data_config_file);
        no_offsets_no_trends.append(stage1obj);

        # The detrended TS without earthquakes or seasonals
        stage2obj = gps_seasonal_removals.make_detrended_ts(stage1obj, 1, 'lssq', data_config_file);
        no_offsets_no_trends_no_seasons.append(stage2obj);

    return [detrended_objects, no_offset_objects, no_offsets_no_trends, no_offsets_no_trends_no_seasons,
            sorted_distances];
