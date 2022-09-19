"""
Viewing a stack of stations
  Step 1: Determine the stations in the radius.
  Step 2: Read them in. Make a list of timeseries objects in one large dataobject.
  Step 3: Compute: Remove outliers, earthquakes, and eventually trend from the data.
  Step 4: Plot in order of increasing latitude, colored by how close they are to the central point
"""

import subprocess
from . import gps_input_pipeline, gps_ts_functions, gps_seasonal_removals, stations_within_radius, offsets
from . import outputs_gps_stacks as out_stack
from file_io import config_io, io_other


def driver(data_config_file, expname, center, radius, proc_center, refframe, outdir, must_include=(None, None)):
    myparams = configure(data_config_file, expname, center, radius, proc_center, refframe, outdir);

    [dataobj_list, offsetobj_list, eqobj_list, paired_distances] = gps_input_pipeline.multi_station_inputs(
        myparams.stations, myparams.blacklist, proc_center, refframe, data_config_file, myparams.distances,
        must_include=must_include);

    [_, no_offset_objects, no_offsets_no_trends, no_offsets_no_trends_no_seasons,
     sorted_distances] = compute(dataobj_list, offsetobj_list, eqobj_list, paired_distances, data_config_file);

    # A series of output options that can be chained together, selected or unselected, etc.
    out_stack.horizontal_full_ts(no_offsets_no_trends, sorted_distances, myparams, "noeq");
    out_stack.horizontal_full_ts(no_offsets_no_trends_no_seasons, sorted_distances, myparams, "noeq_noseasons");
    out_stack.vertical_full_ts(no_offsets_no_trends_no_seasons, sorted_distances, myparams);

    out_stack.horizontal_filtered_plots(no_offsets_no_trends_no_seasons, sorted_distances, myparams);
    out_stack.vertical_filtered_plots(no_offsets_no_trends_no_seasons, sorted_distances, myparams);
    out_stack.vertical_filtered_plots(no_offset_objects, sorted_distances, myparams, 'trendsin_');
    out_stack.pygmt_map(no_offsets_no_trends_no_seasons, myparams);
    return;


def configure(data_config_file, expname, center, radius, proc_center, refframe, outdir):
    # Set up the stacking process
    stations, lons, lats, distances = stations_within_radius.get_stations_within_radius(data_config_file, center,
                                                                                        radius, network=proc_center);
    data_config = config_io.read_config_file(data_config_file);
    blacklist = io_other.read_blacklist(data_config["blacklist"]);
    subprocess.call(["mkdir", "-p", outdir], shell=False);
    outname = expname + "_" + str(center[0]) + "_" + str(center[1]) + "_" + str(radius)
    myparams = out_stack.StackParams(expname=expname, proc_center=proc_center, refframe=refframe, center=center,
                                     radius=radius, stations=stations, bbox=None, distances=distances,
                                     blacklist=blacklist, outdir=outdir, outname=outname, starttime=None, endtime=None,
                                     data_config_file=data_config_file, eqtimes=(), labeltime=None);
    out_stack.write_stack_params(myparams);  # for good measure
    return myparams;


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

    # Objects with no earthquakes or seasonals
    for i in range(len(dataobj_list)):
        # Remove the steps earthquakes
        newobj = offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
        newobj = offsets.remove_offsets(newobj, sorted_eqs[i]);
        newobj = gps_ts_functions.remove_outliers(newobj, 20);  # 20mm outlier definition

        # The detrended TS without earthquakes
        stage1obj = gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq', data_config_file);
        no_offsets_no_trends.append(stage1obj);

        # The detrended TS without earthquakes or seasonals
        stage2obj = gps_seasonal_removals.make_detrended_ts(stage1obj, 1, 'lssq', data_config_file);
        no_offsets_no_trends_no_seasons.append(stage2obj);

    return [detrended_objects, no_offset_objects, no_offsets_no_trends, no_offsets_no_trends_no_seasons,
            sorted_distances];
