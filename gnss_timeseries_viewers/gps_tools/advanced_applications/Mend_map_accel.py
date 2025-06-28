# Python viewing to see the Mendocino stations

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers and earthquakes. Then identify 
# Step 4: Produce a table and plot of accelerations before/after time ranges. 

# Reference:
# Feature: verticals and horizontals at the same time, making two output plots
# Feature: feed seasonal type as parameter, and include that in the output_file name
# This lets us run several experiments.
import gnss_timeseries_viewers.gps_tools.file_io.io_other
import numpy as np
import datetime as dt
import glob, subprocess, os
from gnss_timeseries_viewers.gps_tools import gps_seasonal_removals, offsets, load_gnss
from tectonic_utils.geodesy import haversine
import remove_ets_events

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))  # 5 dirs up
data_config = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt')


def driver(EQcoords, size, network, refframe, fit_type, deltat1, deltat2, expname, station_list=()):
    [stations, outdir, time_after_start_date, critical_variance] = configure(EQcoords, fit_type, size, network,
                                                                             refframe, station_list)
    [dataobj_list, offsetobj_list, eqobj_list] = inputs(stations, network, refframe)
    [noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj] = compute(dataobj_list, offsetobj_list, eqobj_list,
                                                                              deltat1, deltat2, fit_type,
                                                                              time_after_start_date, critical_variance)
    outputs(noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj, outdir, expname, fit_type, network, refframe,
            deltat1, deltat2, time_after_start_date, critical_variance)
    return


def configure(EQcoords, fit_type, overall_size, network, refframe, station_list=()):
    outdir = network + "_" + fit_type + "_" + refframe
    os.makedirs(outdir, exist_ok=True)
    if network == 'nldas' or network == 'gldas' or network == 'noah025' or network == 'lsdm':
        network = 'pbo'  # This is just for finding which stations we will search for.

    time_after_start_date = 7  # optionally, wait a while after the start day.
    critical_variance = 5  # mm/yr. If the time series have a larger variance, we don't consider them

    map_coords = []
    if overall_size == 'medium':
        radius = 350  # km.
    elif overall_size == 'huge':
        radius = -1  # this is a special key for using a coordinate box instead of a radius
        map_coords = [-125.6, -110.0, 32.5, 48.5]
    else:
        radius = 150

    if len(station_list) > 0:
        stations = station_list
    else:
        # Getting the stations of interest ('huge' means we just want within the box.)
        database = load_gnss.create_station_repo(data_config, proc_center=network, refframe=refframe)
        if radius == -1:
            stations = database.search_stations_by_box(map_coords)
        else:

            stations, _ = database.search_stations_by_circle(EQcoords, radius)
        stations = [x.name for x in stations]

    return [stations, outdir, time_after_start_date, critical_variance]


def inputs(station_names, network, refframe):
    database = load_gnss.create_station_repo(data_config, refframe, network)
    dataobj_list, offsetobj_list, eqobj_list = [], [], []
    for station_name in station_names:
        [myData, offset_obj, eq_obj] = database.load_station(station_name)
        if not myData:
            continue
        else:
            dataobj_list.append(myData)
            offsetobj_list.append(offset_obj)
            eqobj_list.append(eq_obj)
    return [dataobj_list, offsetobj_list, eqobj_list]


def add_two_unc_quadrature(unc1, unc2):
    return np.sqrt(unc1 * unc1 + unc2 * unc2)


def compute(dataobj_list, offsetobj_list, eqobj_list, deltat1, deltat2, fit_type, time_after_start, critical_variance):
    dt1_start = dt.datetime.strptime(deltat1[0], "%Y%m%d")
    dt1_end = dt.datetime.strptime(deltat1[1], "%Y%m%d")
    dt2_start = dt.datetime.strptime(deltat2[0], "%Y%m%d")
    dt2_end = dt.datetime.strptime(deltat2[1], "%Y%m%d")

    # No earthquakes objects
    noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj = [], [], [], []

    # For the vertical correction.
    names, coords = [], []
    for i in range(len(dataobj_list)):
        names.append(dataobj_list[i].name)
        coords.append(dataobj_list[i].coords)

    # The main processing loop for slopes.
    for i in range(len(dataobj_list)):
        # Remove the earthquakes
        print(names[i])
        newobj = offsets.remove_offsets(dataobj_list[i], offsetobj_list[i])
        newobj = offsets.remove_offsets(newobj, eqobj_list[i])
        if fit_type == 'none':  # remove seasonal
            newobj = gps_seasonal_removals.make_detrended_ts(newobj, False, fit_type, data_config_file=data_config)
        else:
            if newobj.name == 'P349':
                newobj = gps_seasonal_removals.make_detrended_ts(newobj, True, 'shasta', data_config_file=data_config)
            if newobj.name == 'ORVB':
                newobj = gps_seasonal_removals.make_detrended_ts(newobj, True, 'oroville', data_config_file=data_config)

        if fit_type != 'none':  # remove seasonal
            newobj = gps_seasonal_removals.make_detrended_ts(newobj, True, fit_type, data_config_file=data_config)

        # NOTE: WRITTEN IN JUNE 2019
        # An experiment for removing ETS events
        # if newobj.name in ["P349","P060","P330","P331","P332","P343","P338","P341"]:
        ets_intervals = remove_ets_events.input_tremor_days()
        # newobj=gps_ts_functions.remove_outliers(newobj,3.0)  # 3 mm outlier def.
        # 30 days on either end of the offsets
        # newobj=remove_ets_events.remove_ETS_times(newobj,ets_intervals, offset_num_days=15)
        newobj = remove_ets_events.remove_characteristic_ETS(newobj, ets_intervals)  # using only characteristic offset

        noeq_objects.append(newobj)

        # Get the pre-event and post-event velocities (earthquakes removed)
        starttime = dt1_start + dt.timedelta(days=time_after_start)
        [east_slope_before, north_slope_before, vert_slope_before, esig0, nsig0, _usig0] = newobj.get_slope(starttime,
                                                                                                            dt1_end,
                                                                                                            0.6)
        starttime1 = dt2_start + dt.timedelta(days=time_after_start)
        [east_slope_after, north_slope_after, vert_slope_after, esig1, nsig1, _usig1] = newobj.get_slope(starttime1,
                                                                                                         dt2_end, 0.6)

        # Get the uncertainties on the velocity-change estimate
        starttime2 = dt1_start + dt.timedelta(
            days=time_after_start)
        [east_sl_unc1, north_sl_unc1, vert_sl_unc1] = newobj.get_slope_unc(starttime2, dt1_end)
        starttime3 = dt2_start + dt.timedelta(
            days=time_after_start)
        [east_sl_unc2, north_sl_unc2, vert_sl_unc2] = newobj.get_slope_unc(starttime3, dt2_end)
        east_dv_unc = add_two_unc_quadrature(east_sl_unc1, east_sl_unc2)
        north_dv_unc = add_two_unc_quadrature(north_sl_unc1, north_sl_unc2)
        vert_dv_unc = add_two_unc_quadrature(vert_sl_unc1, vert_sl_unc2)

        # When do we ignore stations? When their detrended time series have a large variance.
        if abs(esig0) > critical_variance or abs(nsig0) > critical_variance or abs(esig1) > critical_variance or abs(
                nsig1) > critical_variance:
            print("Kicking station %s out..." % dataobj_list[i].name)
            [east_slope_after, north_slope_after, vert_slope_after] = [np.nan, np.nan, np.nan]
            [east_slope_before, north_slope_before, vert_slope_before] = [np.nan, np.nan, np.nan]
            [_east_slope_unc1, _north_slope_unc1, _vert_slope_unc1] = [np.nan, np.nan, np.nan]
            [_east_slope_unc2, _north_slope_unc2, _vert_slope_unc2] = [np.nan, np.nan, np.nan]

        east_slope_obj.append([east_slope_before, east_slope_after, east_dv_unc])
        north_slope_obj.append([north_slope_before, north_slope_after, north_dv_unc])
        vert_slope_obj.append([vert_slope_before, vert_slope_after, vert_dv_unc])

    # Adjusting verticals by a reference station.
    vert_slope_obj = vert_adjust_by_reference_stations(names, coords, vert_slope_obj)

    return [noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj]


def vert_adjust_by_reference_stations(names, coords, slope_obj):
    # How do we adjust the verticals for large-scale drought signatures?

    reference_station = 'P208'
    coord_box = [-123, -121, 39, 42]
    eq_coords = [-124.81, 40.53]
    radius = 250
    max_radius = 350
    reference_type = 'radius'  # options = 'radius','box','station'

    new_slope_obj = []
    background_slopes_before = []
    background_slopes_after = []

    for i in range(len(names)):
        if reference_type == 'station':
            if names[i] == reference_station:
                background_slopes_before.append(slope_obj[i][0])
                background_slopes_after.append(slope_obj[i][1])
        elif reference_type == 'box':
            if coord_box[0] < coords[i][0] < coord_box[1]:
                if coord_box[2] < coords[i][1] < coord_box[3]:
                    background_slopes_before.append(slope_obj[i][0])
                    background_slopes_after.append(slope_obj[i][1])
        elif reference_type == 'radius':
            mydistance = haversine.distance([coords[i][1], coords[i][0]], [eq_coords[1], eq_coords[0]])
            if radius < mydistance < max_radius:
                background_slopes_before.append(slope_obj[i][0])
                background_slopes_after.append(slope_obj[i][1])

    vert_reference_before = np.nanmean(background_slopes_before)
    vert_reference_after = np.nanmean(background_slopes_after)
    print("Vert slope before: %f " % vert_reference_before)
    print("Vert slope after: %f " % vert_reference_after)

    for i in range(len(slope_obj)):
        new_slope_obj.append(
            [slope_obj[i][0] - vert_reference_before, slope_obj[i][1] - vert_reference_after, slope_obj[i][2]])

    return new_slope_obj


def outputs(noeq_objects, east_slope_obj, north_slope_obj, vert_slope_obj, outdir, expname, fit_type, network, refframe,
            deltat1, deltat2, time_after_start_date, critical_variance):
    basename = os.path.join(outdir, expname)
    ofile1 = open(basename + '.txt', 'w')
    ofile1.write("# %s network in %s refframe with %s seasonal removal\n" % (network, refframe, fit_type))
    ofile1.write("# %d days gap after EQtime, %s mm/yr maximum variance\n" % (time_after_start_date, critical_variance))
    ofile1.write("# %s minus %s velocities \n" % (deltat2, deltat1))
    for i in range(len(noeq_objects)):
        # Lon, Lat, East, North, 0, Vert, SigE, SigN, SigV, Corr, Name
        ofile1.write("%.2f %.2f %.2f %.2f 0 %.2f %.2f %.2f %.2f 0 %s\n" % (
            noeq_objects[i].coords[0], noeq_objects[i].coords[1], east_slope_obj[i][1] - east_slope_obj[i][0],
            (north_slope_obj[i][1] - north_slope_obj[i][0]), vert_slope_obj[i][1] - vert_slope_obj[i][0],
            east_slope_obj[i][2], north_slope_obj[i][2], vert_slope_obj[i][2], noeq_objects[i].name))
    ofile1.close()

    # Here we call the GMT master script, if we want.
    subprocess.call("./master_plotting.sh " + os.path.join(outdir, expname), shell=True)
    return


# ---------------------- #
# GRACE ONLY FUNCTIONS   # 
# Sometimes we want to see whether loading from GRACE is able to explain the 
# data from Mendocino, etc. 

def grace_driver(deltat1, deltat2, grace_dir, outfile_name, out_dir):
    # For when you want to do the same calculations, for GRACE models
    [file_list, dt1_start, dt1_end, dt2_start, dt2_end, basename] = grace_configure(deltat1, deltat2, grace_dir,
                                                                                    outfile_name)
    [dataobj_list] = grace_inputs(file_list)
    [east_slope_obj, north_slope_obj, vert_slope_obj] = grace_compute(dt1_start, dt1_end, dt2_start, dt2_end,
                                                                      dataobj_list)
    grace_outputs(dataobj_list, east_slope_obj, north_slope_obj, vert_slope_obj, out_dir, basename)
    return


def grace_configure(deltat1, deltat2, grace_dir, outfile_name):
    basename = outfile_name
    dt1_start = dt.datetime.strptime(deltat1[0], "%Y%m%d")
    dt1_end = dt.datetime.strptime(deltat1[1], "%Y%m%d")
    dt2_start = dt.datetime.strptime(deltat2[0], "%Y%m%d")
    dt2_end = dt.datetime.strptime(deltat2[1], "%Y%m%d")

    # Getting the stations of interest
    file_list = glob.glob(grace_dir + "/*.txt")
    print("Reading %d files in %s " % (len(file_list), grace_dir))
    return [file_list, dt1_start, dt1_end, dt2_start, dt2_end, basename]


def grace_inputs(file_list):
    dataobj_list = []
    for item in file_list:
        grace_ts = gnss_timeseries_viewers.gps_tools.file_io.io_other.read_grace(item)
        dataobj_list.append(grace_ts)
    return [dataobj_list]


def grace_compute(dt1_start, dt1_end, dt2_start, dt2_end, dataobject_list):
    east_slope_obj, north_slope_obj, vert_slope_obj = [], [], []
    period_after_start_date = 7  # wait a week.

    for i in range(len(dataobject_list)):
        # Just fit the best line.
        # # Get the pre-event and post-event velocities
        data_ = dataobject_list[i]
        starttime = dt1_start + dt.timedelta(days=period_after_start_date)
        [east_slope_before, north_slope_before, vert_slope_before, _esig0, _nsig0, _usig0] = data_.get_slope(starttime,
                                                                                                             dt1_end,
                                                                                                             0.6)
        data_1 = dataobject_list[i]
        starttime1 = dt2_start + dt.timedelta(days=period_after_start_date)
        [east_slope_after, north_slope_after, vert_slope_after, _esig1, _nsig1, _usig1] = data_1.get_slope(starttime1,
                                                                                                           dt2_end, 0.6)
        east_slope_obj.append([east_slope_before, east_slope_after])
        north_slope_obj.append([north_slope_before, north_slope_after])
        vert_slope_obj.append([vert_slope_before, vert_slope_after])

    # # Experiment: Remove the sinusoidal components. Result is identical.
    # [east_params_before, north_params_before, vert_params_before] = grace_ts_functions.get_linear_annual_semiannual(
    #     dataobject_list[i], starttime=dt1_start+dt.timedelta(days=period_after_start_date),endtime=dt1_end)
    # [east_params_after, north_params_after, vert_params_after] = grace_ts_functions.get_linear_annual_semiannual(
    #     dataobject_list[i],starttime=dt2_start+dt.timedelta(days=period_after_start_date),endtime=dt2_end)
    # east_slope_obj.append([east_params_before[0], east_params_after[0]])
    # north_slope_obj.append([north_params_before[0], north_params_after[0]])
    # vert_slope_obj.append([vert_params_before[0], vert_params_after[0]])

    return [east_slope_obj, north_slope_obj, vert_slope_obj]


def grace_outputs(dataobj_list, east_slope_obj, north_slope_obj, vert_slope_obj, out_dir, basename):
    ofile1 = open(out_dir + basename + '.txt', 'w')
    for i in range(len(dataobj_list)):
        ofile1.write("%f %f %f %f 0 %f 0 0 0 0 %s\n" % (
            dataobj_list[i].coords[0], dataobj_list[i].coords[1], east_slope_obj[i][1] - east_slope_obj[i][0],
            (north_slope_obj[i][1] - north_slope_obj[i][0]), vert_slope_obj[i][1] - vert_slope_obj[i][0],
            dataobj_list[i].name))
    ofile1.close()
    # subprocess.call(['./accel_map_gps.gmt',basename+'.txt',str(map_coords[0]),str(map_coords[1]),str(map_coords[2]),
    #                  str(map_coords[3]),basename],shell=False)
    # print('./accel_map_gps.gmt '+str(map_coords[0])+' '+str(map_coords[1])+' '+str(map_coords[2])
    #       +' '+str(map_coords[3])+' '+basename)
    return
