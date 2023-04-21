# Python viewing of velocities

# Step 1: Determine the stations in the radius. 
# Step 2: Read them in. Make a list of timeseries objects in one large dataobject. 
# Step 3: Compute: Remove outliers, earthquakes, and offsets
# Step 4: Plot maps of vertical velocity


import numpy as np
import collections, subprocess
from GNSS_TimeSeries_Viewers.gps_tools import gps_ts_functions, offsets, load_gnss, vel_functions
import pygmt

Parameters = collections.namedtuple("Parameters",
                                    ['data_config', 'expname', 'proc_center', 'refframe', 'center', 'radius',
                                     'stations', 'outdir', 'outname']);


def driver():
    myparams, database = configure();
    [dataobj_list, offsetobj_list, eqobj_list] = database.load_stations(myparams.stations);
    [_east_slope_list, _north_slope_list, vert_slope_list] = compute(dataobj_list, offsetobj_list, eqobj_list);
    pygmt_vertical_map(myparams, dataobj_list, vert_slope_list);
    return;


def configure():
    # center=[-124, 40.5]; expname='Humboldt'; radius = 100; #
    # center=[-122.5, 40.5]; expname='Chico'; radius = 75; #
    # center=[-119.0, 37.7];     expname='LVC';  radius = 30; # km
    center = [-115.5, 33];
    expname = 'SSGF';
    radius = 80;
    data_config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt";
    blacklist = ["P340", "P316", "P170", "P158", "TRND", "P203", "BBDM", "KBRC", "RYAN", "BEAT", "CAEC",
                 "MEXI"];  # This is global, just keeps growing

    proc_center = 'unr';  # WHICH DATASTREAM DO YOU WANT?
    refframe = 'NA';  # WHICH REFERENCE FRAME?

    database = load_gnss.create_station_repo(data_config_file, proc_center=proc_center, refframe=refframe);
    stations, _ = database.search_stations_by_circle(center, radius);
    stations, _ = vel_functions.remove_blacklist_vels(stations, blacklist);

    outdir = expname + "_" + proc_center
    subprocess.call(["mkdir", "-p", outdir], shell=False);
    outname = expname + "_" + str(center[0]) + "_" + str(center[1]) + "_" + str(radius)
    myparams = Parameters(expname=expname, proc_center=proc_center, refframe=refframe, center=center, radius=radius,
                          stations=[x.name for x in stations], outdir=outdir, outname=outname,
                          data_config=data_config_file);
    return myparams, database;


def compute(dataobj_list, offsetobj_list, eqobj_list):
    east_slope_list, north_slope_list, vert_slope_list = [], [], [];

    # Objects with no earthquakes or seasonals
    for i in range(len(dataobj_list)):
        # Remove the steps earthquakes
        newobj = offsets.remove_offsets(dataobj_list[i], offsetobj_list[i]);
        newobj = offsets.remove_offsets(newobj, eqobj_list[i]);
        [east_slope, north_slope, vert_slope, _east_std, _north_std, _vert_std] = newobj.get_slope(None, None, 0.6);
        east_slope_list.append(east_slope)
        north_slope_list.append(north_slope)
        vert_slope_list.append(vert_slope)

    return [east_slope_list, north_slope_list, vert_slope_list];


def pygmt_vertical_map(myparams, ts_objects, vert_slopes):
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
    min_vert = np.min(vert_slopes);
    if min_vert < -8:
        min_vert = -8;
    max_vert = np.max(vert_slopes);
    if max_vert - min_vert < 5.0:
        label_interval = 0.5;
    else:
        label_interval = 2.0;
    pygmt.makecpt(C="jet", T=str(min_vert - 0.1) + "/" + str(max_vert + 0.1) + "/0.1", H="mycpt.cpt");

    fig = pygmt.Figure()
    fig.basemap(region=region, projection="M8i", B="0.25");
    fig.coast(shorelines="0.5p,black", G='peachpuff2', S='skyblue', D="h");
    fig.coast(N='1', W='1.0p,black');
    fig.coast(N='2', W='0.5p,black');
    fig.text(x=[i + 0.035 for i in lons], y=lats, text=names, font='15p,Helvetica-Bold,black');
    fig.plot(x=geothermals_x, y=geothermals_y, S='i0.2i', G="purple", W='0.5p,black');
    fig.plot(x=lons, y=lats, S='c0.2i', C="mycpt.cpt", G=vert_slopes, W='0.5p,black');
    fig.plot(x=myparams.center[0], y=myparams.center[1], S='a0.1i', G='red', W='0.5p,red')
    fig.colorbar(D="JBC+w4.0i+h", C="mycpt.cpt", G=str(min_vert) + "/" + str(max_vert),
                 B=["x" + str(label_interval), "y+LVelocity(mm/yr)"])
    fig.savefig(myparams.outdir + "/" + myparams.outname + '_map.png');
    return;
