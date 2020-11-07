# Functions for reading a variety of GNSS data formats, including:
# GNSS time series
# GNSS velocity fields
# GRACE and other hydrological models


import numpy as np
import collections, sys, os
import datetime as dt
import configparser

Station_Vel = collections.namedtuple("Station_Vel", ['name', 'nlat', 'elon', 'n', 'e', 'u', 'sn', 'se', 'su',
                                                    'first_epoch',
                                                    'last_epoch', 'refframe', 'proccenter', 'subnetwork', 'survey']);  # in mm/yr, with -180<lon<180, used for velfields
Timeseries = collections.namedtuple("Timeseries", ['name', 'coords', 'dtarray', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su',
                                                   'EQtimes']);  # in mm
Params = collections.namedtuple("Params", ['general_gps_dir', 'pbo_gps_dir', 'unr_gps_dir', 'usgs_gps_dir',
                                           'pbo_earthquakes_dir', 'pbo_offsets_dir',
                                           'unr_offsets_dir', 'unr_coords_file',
                                           'pbo_velocities',
                                           'unr_velocities', 'usgs_vel_dir', 'usgs_networks', 'usgs_cache_file',
                                           'gldas_dir', 'nldas_dir', 'noah_dir', 'grace_dir',
                                           'lsdm_dir', 'stl_dir', 'blacklist']);


def read_config_file(infile):
    if not os.path.isfile(infile):
        print("Error! Data Config file %s not found on your machine. Must fix!" % infile);
        sys.exit(1);

    # Read the all important config file.
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(infile);

    # Where you place all the directories
    general_gps_dir = config.get('py-config', 'gps_data_dir');
    pbo_gps_dir = config.get('py-config', 'pbo_gps_dir');
    unr_gps_dir = config.get('py-config', 'unr_gps_dir');
    usgs_gps_dir = config.get('py-config', 'usgs_gps_dir');
    pbo_earthquakes_dir = config.get('py-config', 'pbo_earthquakes_dir');
    unr_offsets_dir = config.get('py-config', 'unr_offsets_dir');
    pbo_offsets_dir = config.get('py-config', 'pbo_offsets_dir');
    unr_coords_file = config.get('py-config', 'unr_coords_file');
    pbo_velocities = config.get('py-config', 'pbo_velocities');
    unr_velocities = config.get('py-config', 'unr_velocities');
    usgs_vel_dir = config.get('py-config', 'usgs_vel_dir');
    usgs_networks = config.get('py-config', 'usgs_network_list');
    usgs_cache_file = config.get('py-config', 'usgs_cache_file');
    blacklist = config.get('py-config', 'blacklist');
    gldas_dir = config.get('py-config', 'gldas_dir');
    nldas_dir = config.get('py-config', 'nldas_dir');
    noah_dir = config.get('py-config', 'noah_dir');
    grace_dir = config.get('py-config', 'grace_dir');
    lsdm_dir = config.get('py-config', 'lsdm_dir');
    stl_dir = config.get('py-config', 'stl_dir');

    myParams = Params(general_gps_dir=general_gps_dir, pbo_gps_dir=pbo_gps_dir, unr_gps_dir=unr_gps_dir,
                      usgs_gps_dir=usgs_gps_dir, pbo_earthquakes_dir=pbo_earthquakes_dir,
                      pbo_offsets_dir=pbo_offsets_dir, unr_offsets_dir=unr_offsets_dir, unr_coords_file=unr_coords_file,
                      pbo_velocities=pbo_velocities, unr_velocities=unr_velocities, usgs_vel_dir=usgs_vel_dir,
                      usgs_networks=usgs_networks, usgs_cache_file=usgs_cache_file,
                      gldas_dir=gldas_dir, nldas_dir=nldas_dir, noah_dir=noah_dir,
                      grace_dir=grace_dir, lsdm_dir=lsdm_dir, stl_dir=stl_dir, blacklist=blacklist);
    return myParams;


def read_pbo_vel_file(infile):
    # Meant for reading velocity files from the PBO/UNAVCO website.
    # Returns a list of station_vel objects.
    myVelfield = [];
    print("Reading %s" % infile);
    start = 0;
    ifile = open(infile, 'r');
    for line in ifile:
        if start == 1:
            temp = line.split();
            name = temp[0];
            nlat = float(temp[7]);
            elon = float(temp[8]);
            if elon > 180:
                elon = elon - 360.0;
            n = float(temp[19]) * 1000.0;
            e = float(temp[20]) * 1000.0;
            u = float(temp[21]) * 1000.0;
            sn = float(temp[22]) * 1000.0;
            se = float(temp[23]) * 1000.0;
            su = float(temp[24]) * 1000.0;
            t1 = temp[-2];
            t2 = temp[-1];
            first_epoch = dt.datetime.strptime(t1[0:8], '%Y%m%d');
            last_epoch = dt.datetime.strptime(t2[0:8], '%Y%m%d');
            myStationVel = Station_Vel(name=name, elon=elon, nlat=nlat, e=e, n=n, u=u,
                                       se=se, sn=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch,
                                       proccenter='pbo', subnetwork='', refframe='ITRF', survey=0);
            myVelfield.append(myStationVel);
        if "*" in line:
            start = 1;

    ifile.close();
    return [myVelfield];


def read_unr_vel_file(infile, coordinate_file):
    # Meant for reading velocity files from the MAGNET/MIDAS website.
    # Returns a list of station_vel objects.
    myVelfield = []
    # open cache of coordinates and start-dates
    [cache_names, cache_lat, cache_lon, dtbeg, dtend] = np.loadtxt(coordinate_file, unpack=True, skiprows=2,
                                                                   usecols=(0, 1, 2, 7, 8),
                                                dtype={'names': ('name', 'lat', 'lon', 'dtbeg', 'dtend'),
                                                       'formats': ('U4', np.float, np.float, 'U10', 'U10')});

    print("Reading %s" % infile);
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "#":
            continue;
        else:
            name = temp[0];
            e = float(temp[8]) * 1000.0;
            n = float(temp[9]) * 1000.0;
            u = float(temp[10]) * 1000.0;
            se = float(temp[11]) * 1000.0;
            sn = float(temp[12]) * 1000.0;
            su = float(temp[13]) * 1000.0;

            myindex = np.where(cache_names == name);
            if myindex[0] == []:
                print("Error! Could not find cache record for station %s in %s " % (name, coordinate_file));
                sys.exit(0);
            else:
                first_epoch = dt.datetime.strptime(dtbeg[myindex[0][0]], '%Y-%m-%d');
                last_epoch = dt.datetime.strptime(dtend[myindex[0][0]], '%Y-%m-%d');
                elon = cache_lon[myindex[0][0]];
                if elon > 180:
                    elon = elon - 360;
                nlat = cache_lat[myindex[0][0]];

            myStationVel = Station_Vel(name=name, elon=elon, nlat=nlat, e=e, n=n, u=u,
                                       se=se, sn=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch,
                                       proccenter='unr', subnetwork='', refframe='ITRF', survey=0);
            myVelfield.append(myStationVel);
    ifile.close();
    return [myVelfield];


def read_gamit_velfile(infile):
    # Meant for reading a velocity file for example from GAMIT processing
    # Doesn't have starttime and stoptime information.
    myVelfield = [];
    print("Reading %s" % infile);
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "#" or temp[0][0] == "#":
            continue;
        else:
            elon = float(temp[0]);
            if elon > 180:
                elon = elon - 360.0;
            nlat = float(temp[1]);
            e = float(temp[2]);
            n = float(temp[3]);
            se = float(temp[6]);
            sn = float(temp[7]);
            u = float(temp[9]);
            su = float(temp[11]);
            name = temp[12][0:4];
            first_epoch = dt.datetime.strptime("19900101", "%Y%m%d");  # placeholders
            last_epoch = dt.datetime.strptime("20300101", "%Y%m%d");  # placeholders

            myStationVel = Station_Vel(name=name, elon=elon, nlat=nlat, e=e, n=n, u=u,
                               se=se, sn=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch,
                               proccenter='gamit', subnetwork='', refframe='ITRF', survey=0);
            myVelfield.append(myStationVel);
    ifile.close();
    return [myVelfield];


def usgs_vel_file_from_tsfile(infile):
    # Network parsing and plumbing
    # We assume a parallel directory above the TS directory with a bunch of Velocity files
    # From which we can derive coordinates
    usgs_network = infile.split('/')[-2];
    usgs_directory = '';
    for i in range(len(infile.split('/')) - 3):
        usgs_directory = usgs_directory + infile.split('/')[i];
        usgs_directory = usgs_directory + '/';
    network_vel_file = usgs_directory + 'Velocities/NAM_' + usgs_network + '_vels.txt';
    return network_vel_file;


def usgs_network_from_velfile(velfile):
    # velfile is something like 'ITRF_Pacific_Northwest_vels.txt'
    usgs_network = velfile.split('/')[-1][0:-9];
    if usgs_network[0:4] == 'NAM_':
        usgs_network = usgs_network[4:];
        usgs_refframe = 'NAM';
    elif usgs_network[0:5] == 'ITRF_':
        usgs_network = usgs_network[5:];
        usgs_refframe = 'ITRF';
    else:
        usgs_network, usgs_refframe = '', '';
    return usgs_network, usgs_refframe;


def read_usgs_velfile(infile, cache_file):
    # Reading a USGS velocity file.
    # Requires a cache of start and end dates
    print("Reading %s" % infile);
    usgs_sub_network, usgs_refframe = usgs_network_from_velfile(infile);
    [names, lon, lat, e, n, se, sn, u, su] = np.loadtxt(infile, skiprows=3, usecols=(0, 1, 2, 4, 5, 6, 7, 9, 10),
                                                        unpack=True, dtype={'names': (
            'name', 'lon', 'lat', 'evel', 'nvel', 'se', 'sn', 'u', 'su'),
            'formats': ('U4', np.float, np.float, np.float, np.float, np.float,
                np.float, np.float, np.float)})
    # Populating the first_epoch and last_epoch with information from the associated cache.
    [cache_names, startdate, enddate, subnetwork] = np.loadtxt(cache_file, unpack=True, usecols=(0, 3, 4, 5),
                                                               dtype={'names': (
            'name', 'startdate', 'enddate', 'subnetwork'),
            'formats': ('U4', 'U8', 'U8', 'U25')});
    myVelField = [];
    for i in range(len(names)):
        station_name = names[i];
        first_epoch, last_epoch = [None, None];
        idx = np.where(cache_names == station_name);
        num_networks = len(idx[0]);
        for j in range(num_networks):
            single_network = subnetwork[idx[0][j]]
            if single_network == usgs_sub_network:
                first_epoch = dt.datetime.strptime(startdate[idx[0][j]], "%Y%m%d");
                last_epoch = dt.datetime.strptime(enddate[idx[0][j]], "%Y%m%d");
                if first_epoch is None:
                    print("ERROR! Station %s has no first_epoch! Stopping." % station_name ); sys.exit(0);
                if last_epoch is None:
                    print("ERROR! Station %s has no last_epoch! Stopping." % station_name ); sys.exit(0);

        myStationVel = Station_Vel(name=names[i], elon=lon[i], nlat=lat[i], e=e[i], n=n[i], u=u[i],
                                   se=se[i], sn=sn[i], su=su[i], first_epoch=first_epoch, last_epoch=last_epoch,
                                   proccenter='usgs', subnetwork=usgs_sub_network, refframe=usgs_refframe, survey=0);
        myVelField.append(myStationVel);
    return [myVelField];


def read_pbo_pos_file(filename):
    print("Reading %s" % filename);
    [yyyymmdd, Nlat, Elong, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, skiprows=37, unpack=True,
                                                                 usecols=(0, 12, 13, 15, 16, 17, 18, 19, 20));
    dN = [i * 1000.0 for i in dN];
    dE = [i * 1000.0 for i in dE];
    dU = [i * 1000.0 for i in dU];
    Sn = [i * 1000.0 for i in Sn];
    Se = [i * 1000.0 for i in Se];
    Su = [i * 1000.0 for i in Su];
    specific_file = filename.split('/')[-1];
    dtarray = [dt.datetime.strptime(str(int(i)), "%Y%m%d") for i in yyyymmdd];
    myData = Timeseries(name=specific_file[0:4], coords=[Elong[0] - 360, Nlat[0]], dtarray=dtarray, dN=dN, dE=dE, dU=dU,
                        Sn=Sn, Se=Se, Su=Su, EQtimes=[]);
    print("Reading data for station %s in time range %s:%s" % (
        myData.name, dt.datetime.strftime(myData.dtarray[0], "%Y-%m-%d"),
        dt.datetime.strftime(myData.dtarray[-1], "%Y-%m-%d")));
    return [myData];


def read_UNR_magnet_ts_file(filename, coordinates_file):
    print("Reading %s" % filename);
    [datestrs, dE, dN, dU, Se, Sn, Su] = np.loadtxt(filename, usecols=(1, 8, 10, 12, 14, 15, 16), skiprows=1,
                                                    unpack=True,
                                                    dtype={'names': ('datestrs', 'dE', 'dN', 'dU', 'Se', 'Sn', 'Su'),
                                                           'formats': ('U7', np.float, np.float, np.float,
                                                                       np.float, np.float, np.float)});

    station_name = filename.split('/')[-1][0:4];  # the first four characters of the filename
    dtarray = [dt.datetime.strptime(x, '%y%b%d') for x in datestrs];  # has format 07SEP19
    dN = [i * 1000.0 for i in dN];
    dE = [i * 1000.0 for i in dE];
    dU = [i * 1000.0 for i in dU];
    Sn = [i * 1000.0 for i in Sn];
    Se = [i * 1000.0 for i in Se];
    Su = [i * 1000.0 for i in Su];

    [lon, lat] = get_coordinates_for_unr_stations([station_name], coordinates_file);  # format [lat, lon]
    if lon[0] < -360:
        coords = [lon[0] - 360, lat[0]];
    elif lon[0] > 180:
        coords = [lon[0] + 360, lat[0]];
    else:
        coords = [lon[0], lat[0]];

    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    print("Reading data for station %s in time range %s:%s" % (
        myData.name, dt.datetime.strftime(myData.dtarray[0], "%Y-%m-%d"),
        dt.datetime.strftime(myData.dtarray[-1], "%Y-%m-%d")));
    return [myData];


def read_USGS_ts_file(filename):
    print("Reading %s" % filename);
    station_name = filename.split('/')[-1][0:4].upper();  # the first four characters of the filename
    [datestrs, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, unpack=True, usecols=(0, 2, 3, 4, 6, 7, 8),
                                                    dtype={'names': ('datestrs', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su'),
                                                           'formats': (
                                                               'U8', np.float, np.float, np.float, np.float, np.float,
                                                               np.float)})
    dtarray = [dt.datetime.strptime(x, "%Y%m%d") for x in datestrs];
    vel_file = usgs_vel_file_from_tsfile(filename);  # get the associated velocity file for that USGS sub-network
    [names, lon, lat] = np.loadtxt(vel_file, skiprows=3, usecols=(0, 1, 2),
                                                        unpack=True, dtype={'names': (
            'name', 'lon', 'lat'), 'formats': ('U4', np.float, np.float)})
    coords = None;   # default value in case station wasn't found
    for i in range(len(names)):
        if names[i] == station_name:
            coords = [lon[i], lat[i]];
    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    return [myData];


def read_pbo_hydro_file(filename, coords_file=None):
    # Useful for reading hydrology files like NLDAS, GLDAS, etc.
    # In the normal pipeline for this function, it is guaranteed to be given a real file.
    print("Reading %s" % filename);
    dtarray = [];
    station_name = filename.split('/')[-1][0:4];
    station_name = station_name.upper();
    if coords_file is not None:
        [lon, lat] = get_coordinates_for_unr_stations([station_name], coords_file);  # format [lat, lon]
        coords = [lon[0], lat[0]];
    else:
        coords = [None, None];  # can return an object without meaningful coordinates if not asked for them.
    [dts, dN, dE, dU] = np.loadtxt(filename, usecols=(0, 3, 4, 5), dtype={'names': ('dts', 'dN', 'dE', 'dU'),
                                                                          'formats': (
                                                                              'U10', np.float, np.float, np.float)},
                                   skiprows=20, delimiter=',', unpack=True);
    for i in range(len(dts)):
        dtarray.append(dt.datetime.strptime(dts[i], "%Y-%m-%d"));
    Se = 0.2 * np.ones(np.shape(dE));
    Sn = 0.2 * np.ones(np.shape(dN));
    Su = 0.2 * np.ones(np.shape(dU));

    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    return [myData];


def read_lsdm_file(filename, coords_file=None):
    # Useful for reading hydrology files from LSDM German loading product
    # In the normal pipeline for this function, it is guaranteed to be given a real file.
    print("Reading %s" % filename);
    dtarray = [];
    station_name = filename.split('/')[-1][0:4];
    if coords_file is not None:
        [lon, lat] = get_coordinates_for_unr_stations([station_name], coords_file);  # format [lat, lon]
        coords = [lon[0], lat[0]];
    else:
        coords = [None, None];  # can return an object without meaningful coordinates if not asked for them.
    [dts, dU, dN, dE] = np.loadtxt(filename, usecols=(0, 1, 2, 3), dtype={'names': ('dts', 'dN', 'dE', 'dU'),
                                                                          'formats': (
                                                                              'U10', np.float, np.float, np.float)},
                                   skiprows=3, delimiter=',', unpack=True);
    for i in range(len(dts)):
        dtarray.append(dt.datetime.strptime(dts[i][0:10], "%Y-%m-%d"));
    dN = [i * 1000.0 for i in dN];
    dE = [i * 1000.0 for i in dE];
    dU = [i * 1000.0 for i in dU];
    Se = 0.2 * np.ones(np.shape(dE));
    Sn = 0.2 * np.ones(np.shape(dN));
    Su = 0.2 * np.ones(np.shape(dU));
    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    return [myData];


def get_coordinates_for_unr_stations(station_names, coordinates_file):
    # station_names is an array
    # open cache of coordinates and start-dates
    [cache_names, cache_lat, cache_lon] = np.loadtxt(coordinates_file, unpack=True, skiprows=2,
                                                                   usecols=(0, 1, 2),
                                                dtype={'names': ('name', 'lat', 'lon'),
                                                       'formats': ('U4', np.float, np.float)});
    lon, lat = [], [];
    # find the stations
    for i in range(len(station_names)):
        myindex = cache_names.index(station_names[i]);
        if myindex == []:
            print("Error! Could not find UNR coordinates for station %s " % station_names[i]);
            print("Returning [0,0]. ");
            lon.append(0.0);
            lat.append(0.0);
        else:
            elon = cache_lon[myindex];
            if elon > 180:
                elon = elon - 360;
            lon.append(elon);
            lat.append(cache_lat[myindex]);
    return [lon, lat];


def read_grace(filename):
    # Read the GRACE data into a GPS-style time series object.
    # THE GRACE DATA
    print("Reading %s" % filename);
    station_name = filename.split('/')[-1];  # this is the local name of the file
    station_name = station_name.split('_')[1];  # this is the four-character name
    try:
        [dts, lon, lat, _, u, v, w] = np.loadtxt(filename, usecols=range(0, 7),
                                                 dtype={'names': ('dts', 'lon', 'lat', 'temp', 'u', 'v', 'w'),
                                                        'formats': (
                                                            'U11', np.float, np.float, np.float, np.float, np.float,
                                                            np.float)}, unpack=True);
    except FileNotFoundError:
        print("ERROR! Cannot find GRACE model for file %s" % filename);
        return None;
    grace_t = [];
    for i in range(len(dts)):
        grace_t.append(dt.datetime.strptime(dts[i], "%d-%b-%Y"));
    grace_t = [i + dt.timedelta(days=15) for i in
               grace_t];  # we add 15 days to plot the GRACE data at the center of the bin.
    u = np.array(u);
    v = np.array(v);
    w = np.array(w);
    S = np.zeros(np.shape(u));
    GRACE_TS = Timeseries(name=station_name, coords=[lon[0], lat[0]], dtarray=grace_t, dE=u, dN=v, dU=w, Se=S, Sn=S,
                          Su=S, EQtimes=[]);
    return [GRACE_TS];


# ---------- WRITING FUNCTIONS ----------- # 
def write_pbo_pos_file(ts_object, filename, comment=""):
    # Useful for writing common mode objects, etc.
    # Opposite of read_pbo_pos_file(filename)
    print("Writing %s" % filename);
    ofile = open(filename, 'w');
    ofile.write("%s\n" % comment);
    for i in range(36):
        ofile.write("/\n");
    for i in range(len(ts_object.dtarray)):
        ofile.write("%s 0 0 0 0 0 0 0 0 0 0 0 " % (
            dt.datetime.strftime(ts_object.dtarray[i], "%Y%m%d")));  # the first 12 columns
        ofile.write("%.5f %.5f 0 " % (ts_object.coords[1], ts_object.coords[0] + 360));
        ofile.write("%.6f %.6f %.6f %.6f %.6f %.6f\n" % (
            ts_object.dN[i] / 1000, ts_object.dE[i] / 1000, ts_object.dU[i] / 1000.0, ts_object.Sn[i] / 1000.0,
            ts_object.Se[i] / 1000.0, ts_object.Su[i] / 1000.0));
    ofile.close();
    return;


def write_humanread_vel_file(myVelfield, outfile):
    ofile = open(outfile, 'w');
    ofile.write(
        "Format: lon(deg) lat(deg) e(mm) n(mm) u(mm) Se(mm) Sn(mm) Su(mm) first_dt(yyyymmdd) last_dt(yyyymmdd) name\n");
    for station_vel in myVelfield:
        first_epoch = dt.datetime.strftime(station_vel.first_epoch, '%Y%m%d');
        last_epoch = dt.datetime.strftime(station_vel.last_epoch, '%Y%m%d');
        ofile.write("%f %f %f %f %f %f %f %f %s %s %s\n" % (
            station_vel.elon, station_vel.nlat, station_vel.e, station_vel.n, station_vel.u, station_vel.se,
            station_vel.sn, station_vel.su, first_epoch, last_epoch, station_vel.name));
    ofile.close();
    return;


def write_gmt_velfile(myVelfield, outfile):
    ofile = open(outfile, 'w');
    ofile.write("# Format: lon(deg) lat(deg) e(mm) n(mm) Se(mm) Sn(mm) 0 0 1 name\n");
    for station_vel in myVelfield:
        if station_vel.sn < 0.2:  # trying to make a clean dataset
            ofile.write("%f %f %f %f %f %f 0 0 1 %s\n" % (
                station_vel.elon, station_vel.nlat, station_vel.e, station_vel.n, station_vel.se,
                station_vel.sn, station_vel.name));
    ofile.close();
    return;
