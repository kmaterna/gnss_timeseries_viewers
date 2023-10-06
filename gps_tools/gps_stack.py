"""
Viewing a stack of stations
  Step 1: Determine the stations in the radius.
  Step 2: Read them in. Make a list of timeseries objects in one large dataobject.
  Step 3: Compute: Remove outliers, earthquakes, and trend from the data.
  Step 4: Plot in order of increasing latitude, colored by how close they are to the central point
"""
import os
from . import gps_seasonal_removals, offsets, load_gnss, vel_functions, pygmt_plots, utilities
from . import outputs_gps_stacks as out_stack
from .file_io import config_io, io_other, io_nota
import datetime as dt


def driver(data_config, expname, center, radius, proc_center, refframe, outdir,
           starttime=dt.datetime.strptime("20050505", "%Y%m%d"), endtime=None, write_pos_files=False):
    """
    :param data_config: string, path to data config file
    :param expname: string, used in outdir name
    :param center: (lon, lat)
    :param radius: float, km
    :param proc_center: string
    :param refframe: string
    :param outdir: string
    :param starttime: dt.datetime
    :param endtime: dt.datetime
    :param write_pos_files: boolean to decide whether to write output files
    """
    outname = configure(expname, center, radius, outdir);

    database, stations, distances = build_database(data_config, proc_center, refframe, center, radius);
    [dataobj_list, offsetobj_list, eqobj_list] = database.load_stations(stations);
    param_dict = pack_params(data_config, expname, center, radius, proc_center, refframe, outdir,
                             starttime, endtime, stations);   # for tracking metadata
    out_stack.write_params(outfile=os.path.join(outdir, outname + "_stack_params.txt"), param_dict=param_dict);
    [detr, _, no_offsets_detr, no_offsets_detr_deseas, distances] = compute(dataobj_list, offsetobj_list, eqobj_list,
                                                                            distances, data_config);

    # A series of output options that can be chained together, selected or unselected, etc.
    out_stack.horizontal_full_ts(detr, distances, outname=os.path.join(outdir, outname+'_TS.png'),
                                 start_time_plot=starttime, end_time_plot=endtime);
    out_stack.horizontal_full_ts(no_offsets_detr, distances, outname=os.path.join(outdir, outname+'_TS_noeq.png'),
                                 start_time_plot=starttime, end_time_plot=endtime);
    out_stack.horizontal_full_ts(no_offsets_detr_deseas, distances, start_time_plot=starttime, end_time_plot=endtime,
                                 outname=os.path.join(outdir, outname + '_TS_noeq_noseasons.png'));
    out_stack.vertical_full_ts(no_offsets_detr_deseas, distances, start_time_plot=starttime, end_time_plot=endtime,
                               outname=os.path.join(outdir, outname + '_TS_vertical.png'));

    out_stack.horizontal_filtered_plots(no_offsets_detr_deseas, distances, start_time_plot=starttime,
                                        outname=os.path.join(outdir, outname + '_TS_horiz_filt.png'),
                                        end_time_plot=endtime);
    out_stack.vertical_filtered_plots(no_offsets_detr_deseas, distances,
                                      outname=os.path.join(outdir, outname + '_TS_vert_detrended_filt.png'),
                                      start_time_plot=starttime, end_time_plot=endtime);

    pygmt_plots.map_ts_objects(dataobj_list, os.path.join(outdir, outname+'_map.png'), center=center);
    if write_pos_files:
        write_stack_pos_files(detr, outdir, outname);
    return no_offsets_detr_deseas;


def build_database(data_config_file, proc_center, refframe, center, radius):
    # Set up the stacking process, return list of desired stations by name
    database = load_gnss.create_station_repo(data_config_file, proc_center=proc_center, refframe=refframe);
    station_names, _ = database.search_stations_by_circle(center, radius, basic_clean=True);
    data_config = config_io.read_config_file(data_config_file);
    blacklist = io_other.read_blacklist(data_config["blacklist"]);
    station_names = vel_functions.remove_blacklist_vels(station_names, blacklist);
    return database, [x.name for x in station_names], [x.get_distance_to_point(center) for x in station_names];


def configure(expname, center, radius, outdir):
    os.makedirs(outdir, exist_ok=True);
    outname = expname + "_" + str(center[0]) + "_" + str(center[1]) + "_" + str(radius)
    return outname;


def pack_params(data_config_file, expname, center, radius, proc_center, refframe, outdir, starttime, endtime, stations):
    return locals();


def compute(dataobj_list, offsetobj_list, eqobj_list, distances, data_config_file):
    latitudes_list = [i.coords[1] for i in dataobj_list];
    sorted_objects = [x for _, x in sorted(zip(latitudes_list, dataobj_list))];  # the raw, sorted data.
    sorted_offsets = [x for _, x in sorted(zip(latitudes_list, offsetobj_list))];  # the raw, sorted data.
    sorted_eqs = [x for _, x in sorted(zip(latitudes_list, eqobj_list))];  # the raw, sorted data.
    sorted_distances = [x for _, x in sorted(zip(latitudes_list, distances))];  # the sorted distances.

    detr_objs, no_offset_objs, no_offsets_detr, no_offsets_detr_deseas = [], [], [], [];

    # Detrended objects (or objects with trends and no offsets; depends on what you want.)
    for i in range(len(sorted_objects)):
        newobj = gps_seasonal_removals.make_detrended_ts(sorted_objects[i], 0, 'lssq', data_config_file);
        newobj = newobj.remove_outliers(20);  # 20mm outlier definition
        newobj = offsets.remove_offsets(newobj, sorted_offsets[i]);  # mechanical offsets
        detr_objs.append(newobj);  # still has offsets, doesn't have trends

        newobj = offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
        newobj = offsets.remove_offsets(newobj, sorted_eqs[i]);
        no_offset_objs.append(newobj);  # still has trends, doesn't have offsets

    # # Objects with no earthquakes or seasonals
    for i in range(len(sorted_objects)):
        # Remove the steps earthquakes
        newobj = offsets.remove_offsets(sorted_objects[i], sorted_offsets[i]);
        newobj = offsets.remove_offsets(newobj, sorted_eqs[i]);
        newobj = newobj.remove_outliers(20);  # 20mm outlier definition

        # The detrended TS without earthquakes
        stage1obj = gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq', data_config_file);
        no_offsets_detr.append(stage1obj);

        # The detrended TS without earthquakes or seasonals
        stage2obj = gps_seasonal_removals.make_detrended_ts(stage1obj, 1, 'lssq', data_config_file);
        no_offsets_detr_deseas.append(stage2obj);

    return [detr_objs, no_offset_objs, no_offsets_detr, no_offsets_detr_deseas, sorted_distances];


def extract_stack_displacements(data_obj_list, starttime, endtime, num_days=3):
    """
    Extract displacements between starttime and endtime from a stack of Timeseries objects.
    Return as a list of Station_Vel objects.
    Uncertainties on the displacements are also computed by added in quadrature.
    """
    displacements = [];
    for ts in data_obj_list:
        e1, n1, u1 = ts.get_values_at_date(starttime, num_days=num_days);
        e2, n2, u2 = ts.get_values_at_date(endtime, num_days=num_days);
        se1, sn1, su1 = ts.get_uncertainties_for_avg_position(starttime, num_days);
        se2, sn2, su2 = ts.get_uncertainties_for_avg_position(endtime, num_days);
        newpt = vel_functions.Station_Vel(elon=ts.coords[0], nlat=ts.coords[1], e=(e2-e1), n=(n2-n1), u=(u2-u1),
                                          se=utilities.add_in_quadrature(se1, se2),
                                          sn=utilities.add_in_quadrature(sn1, sn2),
                                          su=utilities.add_in_quadrature(su1, su2),
                                          first_epoch=starttime, last_epoch=endtime, name=ts.name);
        displacements.append(newpt);
    return displacements;


def write_stack_pos_files(data_obj_list, outdir, outname):
    os.makedirs(os.path.join(outdir, 'TS_'+outname), exist_ok=True);
    for item in data_obj_list:
        io_nota.write_pbo_pos_file(item, os.path.join(outdir, 'TS_'+outname, item.name+'_detrended.pos'));
    return;
