# The purpose of this script is to compute a common mode term (average of the time series)
# We save that as a time series format
# We plot the vertical time series with and without the common mode term.

import gnss_timeseries_viewers.gps_tools.file_io.io_nota as io_nota
from gnss_timeseries_viewers.gps_tools.gps_ts_functions import Timeseries
import numpy as np
import scipy.ndimage
import collections, os
import datetime as dt
from gnss_timeseries_viewers.gps_tools import gps_seasonal_removals, \
    offsets, outputs_gps_stacks, load_gnss, vel_functions
import pygmt

Parameters = collections.namedtuple("Parameters",
                                    ['data_config', 'expname', 'proc_center', 'refframe', 'center', 'radius',
                                     'stations', 'distances', 'blacklist', 'outdir', 'outname']);


def driver():
    myparams, database = configure();
    [dataobj_list, offsetobj_list, eqobj_list] = database.load_stations(myparams.stations);
    [common_mode, raw_objects, cmr_objects, deltas] = compute(dataobj_list, offsetobj_list, eqobj_list,
                                                              myparams.data_config);
    print(deltas);
    outputs_gps_stacks.vertical_filtered_plots(raw_objects, (), myparams, "vertical_filt");
    outputs_gps_stacks.vertical_filtered_plots(cmr_objects, (), myparams, "no_cmr_filt");
    pygmt_map(cmr_objects, myparams, deltas);
    write_cm_object(common_mode, myparams);
    return;


def configure():
    # center=[-125.134, 40.829]; expname='Mend'; radius = 120; # Keeper for the 2014 earthquake
    # center=[-125.134, 40.829]; expname='Mend'; radius = 90; #
    # center=[-124, 40.5]; expname='Humboldt'; radius = 60; #
    # center=[-122.5, 40.5]; expname='Chico'; radius = 75; #
    # center=[-124.0, 38.0];     expname='Nbay'; radius = 125;
    # center=[-119.0, 34.5];     expname='SoCal';  radius = 25; # km
    # center=[-116.0, 34.5];     expname='Mojave';  radius = 35; # km
    # center=[-117.5, 35.5];     expname='ECSZ';  radius = 50; # km
    # center=[-119.0, 37.7];     expname='LVC';  radius = 30; # km
    center = [-115.5, 32.85];
    expname = 'SSGF';
    radius = 30;
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
    data_config_file = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt');
    blacklist = ["P316", "P170", "P158", "TRND", "P203", "BBDM", "KBRC", "RYAN", "BEAT", "CAEC", "MEXI", "BOMG",
                 "FSHB"];  # This is global, just keeps growing
    # center=[-115.5, 33]; expname='SSGF'; radius = 40;

    proc_center = 'cwu';  # WHICH DATASTREAM DO YOU WANT?
    refframe = 'NA';  # WHICH REFERENCE FRAME?

    database = load_gnss.create_station_repo(data_config_file, proc_center=proc_center, refframe=refframe);
    stations, _ = database.search_stations_by_circle(center, radius);
    stations, _ = vel_functions.remove_blacklist_vels(stations, blacklist);

    outdir = expname + "_" + proc_center
    os.makedirs(outdir, exist_ok=True);
    outname = expname + "_" + str(center[0]) + "_" + str(center[1]) + "_" + str(radius)
    myparams = Parameters(data_config=data_config_file,
                          expname=expname, proc_center=proc_center, refframe=refframe, center=center, radius=radius,
                          stations=[x.name for x in stations], distances=(),
                          blacklist=blacklist, outdir=outdir, outname=outname);
    return myparams, database;


def compute(dataobj_list, offsetobj_list, eqobj_list, data_config_file):
    print("Computing common mode from %d stations " % (len(dataobj_list)))
    detrended_objects = [];
    cmr_objects = [];  # common-mode-removed objects
    cmr_deviations = [];
    # Objects with no earthquakes or seasonals
    for i in range(len(dataobj_list)):
        newobj = offsets.remove_offsets(dataobj_list[i], offsetobj_list[i]);
        newobj = offsets.remove_offsets(newobj, eqobj_list[i]);
        newobj = gps_seasonal_removals.make_detrended_ts(newobj, 0, 'lssq', data_config_file);
        detrended_objects.append(newobj);

    common_mode_obj = define_common_mode(detrended_objects);

    print("Removing common mode from %d objects " % (len(dataobj_list)))
    for i in range(len(dataobj_list)):
        cmr_object = remove_cm_from_object(detrended_objects[i], common_mode_obj);
        cmr_objects.append(cmr_object);

    # Keep a measure of the spread of this data
    for i in range(len(cmr_objects)):
        udata = scipy.ndimage.median_filter(cmr_objects[i].dU, size=365);
        cmr_deviations.append(np.nanmax(udata) - np.nanmin(udata));

    return [common_mode_obj, detrended_objects, cmr_objects, cmr_deviations];


def define_common_mode(detrended_objects):
    print("Defining common mode from %d objects" % (len(detrended_objects)))
    total_dtarray, data_array_dE, data_array_dN, data_array_dU = make_ND_arrays(detrended_objects)

    cm_dn, cm_de, cm_du, cm_sn, cm_se, cm_su = [], [], [], [], [], [];
    for i in range(len(total_dtarray)):
        cm_dn.append(get_common_mode_from_list(data_array_dN[i, :]))
        cm_de.append(get_common_mode_from_list(data_array_dE[i, :]))
        cm_du.append(get_common_mode_from_list(data_array_dU[i, :]))
    sigmas = np.array([1 for _i in cm_dn]);
    common_mode_obj = Timeseries(name='como', coords=[-115, 32], dtarray=np.array(total_dtarray), dN=np.array(cm_dn),
                                 dE=np.array(cm_de), dU=np.array(cm_du), Se=sigmas, Sn=sigmas, Su=sigmas, EQtimes=[]);
    return common_mode_obj;


def get_common_mode_from_list(data_list):
    # Given a list of stations' values for a certain day, what is the common mode?
    # Taking median, since it's robust to outliers
    if sum(~np.isnan(data_list)) < 5:
        value = 0;
    else:
        value = np.nanmedian(data_list);
    return value;


def remove_cm_from_object(data_object, cm_object):
    cmr_dE, cmr_dN, cmr_dU = [], [], [];
    for i in range(len(data_object.dtarray)):
        if data_object.dtarray[i] in cm_object.dtarray:
            cm_object_indices = [d == data_object.dtarray[i] for d in
                                 cm_object.dtarray];  # select which day you're doing
            cm_value = cm_object.dE[cm_object_indices];  # get the index of the cm_object's day
            cmr_dE.append(data_object.dE[i] - cm_value);
            cm_value = cm_object.dN[cm_object_indices];
            cmr_dN.append(data_object.dN[i] - cm_value);
            cm_value = cm_object.dU[cm_object_indices];
            cmr_dU.append(data_object.dU[i] - cm_value);

    sigmas = np.array([1 for _i in cmr_dU]);
    cmr_object = Timeseries(name=data_object.name, coords=data_object.coords, dtarray=data_object.dtarray,
                            dN=np.array(cmr_dN), dE=np.array(cmr_dE), dU=np.array(cmr_dU),
                            Se=sigmas, Sn=sigmas, Su=sigmas, EQtimes=data_object.EQtimes);
    return cmr_object;


def make_ND_arrays(object_list):
    # Make a 2D array for each of the components.
    first_date = object_list[0].dtarray[0];
    last_date = object_list[0].dtarray[-1];
    total_dtarray = [];
    # Get the list of
    for i in range(len(object_list)):
        if object_list[i].dtarray[0] < first_date:
            first_date = object_list[i].dtarray[0];
        if object_list[i].dtarray[-1] > last_date:
            last_date = object_list[i].dtarray[-1];
    total_dtarray.append(first_date);

    for i in range(1, 10000):
        if total_dtarray[-1] < last_date:
            total_dtarray.append(first_date + dt.timedelta(days=i));

    data_array_dE = np.full([len(total_dtarray), len(object_list)], np.nan);
    data_array_dN = np.full([len(total_dtarray), len(object_list)], np.nan);
    data_array_dU = np.full([len(total_dtarray), len(object_list)], np.nan);

    # Here's a bit slow.
    for i in range(len(total_dtarray)):  # for each day...
        for j in range(len(object_list)):  # check each station
            indices = [d == total_dtarray[i] for d in object_list[j].dtarray];  # select which day you're doing
            # A boolean list.
            if True in indices:
                data_array_dE[i][j] = object_list[j].dE[indices]  # please don't have a duplicated day...
                data_array_dN[i][j] = object_list[j].dN[indices]  # please don't have a duplicated day...
                data_array_dU[i][j] = object_list[j].dU[indices]  # please don't have a duplicated day...

    return total_dtarray, data_array_dE, data_array_dN, data_array_dU;


def pygmt_map(ts_objects, myparams, deltas):
    offset = 0.2;
    geothermals_x = [-115.528300, -115.250000, -115.515300, -115.600000];
    geothermals_y = [32.716700, 32.783300, 33.015300, 33.200000];

    lons, lats, names = [], [], [];
    for i in range(len(ts_objects)):
        lons.append(ts_objects[i].coords[0]);
        lats.append(ts_objects[i].coords[1]);
        names.append(ts_objects[i].name);
    region = [min(lons) - offset, max(lons) + offset, min(lats) - offset, max(lats) + offset];

    # Make a new color bar
    min_vert = np.min(deltas);
    _max_vert = np.max(deltas);
    max_vert = 20;
    label_interval = 2.0;
    pygmt.makecpt(C="jet", T=str(min_vert - 0.1) + "/" + str(max_vert + 0.1) + "/0.1", H="mycpt.cpt");

    fig = pygmt.Figure()
    fig.basemap(region=region, projection="M8i", B="0.25");
    fig.coast(shorelines="0.5p,black", G='peachpuff2', S='skyblue', D="h");
    fig.coast(N='1', W='1.0p,black');
    fig.coast(N='2', W='0.5p,black');
    fig.text(x=[i + 0.035 for i in lons], y=lats, text=names, font='15p,Helvetica-Bold,black');
    fig.plot(x=geothermals_x, y=geothermals_y, S='i0.2i', G="purple", W='0.5p,black');
    fig.plot(x=lons, y=lats, S='c0.2i', C="mycpt.cpt", G=deltas, W='0.5p,black');
    fig.plot(x=myparams.center[0], y=myparams.center[1], S='a0.1i', G='red', W='0.5p,red')
    fig.colorbar(D="JBC+w4.0i+h", C="mycpt.cpt", G=str(min_vert) + "/" + str(max_vert),
                 B=["x" + str(label_interval), "y+LVert(mm)"])
    fig.savefig(myparams.outdir + "/" + myparams.outname + '_map.png');
    return;


def write_cm_object(common_mode, myparams):
    filename = os.path.join(myparams.outdir, myparams.outname + "_common_mode.pos");
    print("Writing common_mode into %s " % filename)
    comment = myparams.stations;  # right now this is all stations, not just the ones that were used.
    io_nota.write_pbo_pos_file(common_mode, filename, comment);
    return;
