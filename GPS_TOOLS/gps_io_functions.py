# Functions for reading a variety of GNSS data formats, including:
# GNSS time series
# GNSS velocity fields
# GRACE and other hydrological models


import numpy as np
import collections, sys, os
import datetime as dt
import configparser

Velfield = collections.namedtuple("Velfield", ['name', 'nlat', 'elon', 'n', 'e', 'u', 'sn', 'se', 'su', 'first_epoch',
                                               'last_epoch']);  # in mm/yr, with -180<lon<180
Timeseries = collections.namedtuple("Timeseries", ['name', 'coords', 'dtarray', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su',
                                                   'EQtimes']);  # in mm
Params = collections.namedtuple("Params", ['general_gps_dir','pbo_gps_dir', 'unr_gps_dir', 'usgs_gps_dir', 
                                           'pbo_earthquakes_dir', 'pbo_offsets_dir',
                                           'unr_offsets_dir', 'unr_coords_file', 'pbo_velocities',
                                           'unr_velocities', 'usgs_velocities', 'usgs_networks',
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
    usgs_velocities = config.get('py-config', 'usgs_vel_dir');
    usgs_networks = config.get('py-config', 'usgs_network_list');
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
                      pbo_velocities=pbo_velocities, unr_velocities=unr_velocities, usgs_velocities=usgs_velocities,
                      usgs_networks=usgs_networks,
                      gldas_dir=gldas_dir, nldas_dir=nldas_dir, noah_dir=noah_dir,
                      grace_dir=grace_dir, lsdm_dir=lsdm_dir, stl_dir=stl_dir, blacklist=blacklist);
    return myParams;


def read_pbo_vel_file(infile):
    # Meant for reading velocity files from the PBO/UNAVCO website.
    # Returns a Velfield object.
    print("Reading %s" % infile);
    start = 0;
    ifile = open(infile, 'r');
    name = [];
    nlat, elon = [], [];
    n, e, u = [], [], [];
    sn, se, su = [], [], [];
    first_epoch = [];
    last_epoch = [];
    for line in ifile:
        if start == 1:
            temp = line.split();
            name.append(temp[0]);
            nlat.append(float(temp[7]));
            elon_temp = float(temp[8]);
            if elon_temp > 180:
                elon_temp = elon_temp - 360.0;
            elon.append(elon_temp);
            n.append(float(temp[19]) * 1000.0);
            e.append(float(temp[20]) * 1000.0);
            u.append(float(temp[21]) * 1000.0);
            sn.append(float(temp[22]) * 1000.0);
            se.append(float(temp[23]) * 1000.0);
            su.append(float(temp[24]) * 1000.0);
            t1 = temp[-2];
            t2 = temp[-1];
            first_epoch.append(dt.datetime.strptime(t1[0:8], '%Y%m%d'));
            last_epoch.append(dt.datetime.strptime(t2[0:8], '%Y%m%d'));
        if "*" in line:
            start = 1;
    ifile.close();
    myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch,
                          last_epoch=last_epoch);
    return [myVelfield];


def read_unr_vel_file(infile, coordinate_file):
    # Meant for reading velocity files from the MAGNET/MIDAS website.
    # Returns a Velfield object.
    print("Reading %s" % infile);
    name = [];
    n, e, u = [], [], [];
    sn, se, su = [], [], [];
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "#":
            continue;
        else:
            name.append(temp[0]);
            e.append(float(temp[8]) * 1000.0);
            n.append(float(temp[9]) * 1000.0);
            u.append(float(temp[10]) * 1000.0);
            se.append(float(temp[11]) * 1000.0);
            sn.append(float(temp[12]) * 1000.0);
            su.append(float(temp[13]) * 1000.0);
    ifile.close();

    [elon, nlat] = get_coordinates_for_stations(name, coordinate_file);
    [first_epoch, last_epoch] = get_start_times_for_unr_stations(name, coordinate_file);

    myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch,
                          last_epoch=last_epoch);
    return [myVelfield];


def read_gamit_velfile(infile):
    # Meant for reading a velocity file for example from GAMIT processing
    # Doesn't have starttime and stoptime information.
    print("Reading %s" % infile);
    ifile = open(infile, 'r');
    name = [];
    nlat, elon = [], [];
    n, e, u = [], [], [];
    sn, se, su = [], [], [];
    first_epoch = [];
    last_epoch = [];

    for line in ifile:
        temp = line.split();
        if temp[0] == "#" or temp[0][0] == "#":
            continue;
        else:
            elon_temp = float(temp[0]);
            if elon_temp > 180:
                elon_temp = elon_temp - 360.0;
            elon.append(elon_temp);
            nlat.append(float(temp[1]));
            e.append(float(temp[2]));
            n.append(float(temp[3]));
            se.append(float(temp[6]));
            sn.append(float(temp[7]));
            u.append(float(temp[9]));
            su.append(float(temp[11]));
            name.append(temp[12][0:4]);
            first_epoch.append(dt.datetime.strptime("19900101", "%Y%m%d"));  # placeholders
            last_epoch.append(dt.datetime.strptime("20300101", "%Y%m%d"));  # placeholders

    myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch,
                          last_epoch=last_epoch);

    return [myVelfield];


def read_usgs_velfile(infile, ts_directory=None):
    # Reading a USGS velocity file.
    # If we provide a ts_directory, we can get first_epoch and last_epoch for each station.
    # Otherwise, we have to put placeholders.
    print("Reading %s" % infile);
    [name, lon, lat, e, n, se, sn, u, su] = np.loadtxt(infile, skiprows=3, usecols=(0, 1, 2, 4, 5, 6, 7, 9, 10),
                                                       unpack=True, dtype={'names': (
            'name', 'lon', 'lat', 'evel', 'nvel', 'se', 'sn', 'u', 'su'),
            'formats': (
                'U4', np.float, np.float, np.float, np.float, np.float,
                np.float, np.float, np.float)})
    # Populating the first_epoch and last_epoch with information from the associated time series directory.
    first_epoch = [];
    last_epoch = [];
    if ts_directory is not None:
        for station in name:
            ts_filename = ts_directory + '/' + station + '_NAfixed.rneu';
            [dates, _] = np.loadtxt(ts_filename, unpack=True, usecols=(0, 1),
                                    dtype={'names': ('dtstrs', 'datenums'), 'formats': ('U8', np.float)});
            if np.size(dates) == 1:  # if the campaign station only had one day of data...
                first_epoch = [dt.datetime.strptime(str(dates), "%Y%m%d") for x in name];
                last_epoch = [dt.datetime.strptime(str(dates), "%Y%m%d") for x in name];
            else:
                first_epoch.append(dt.datetime.strptime(dates[0], "%Y%m%d"));
                last_epoch.append(dt.datetime.strptime(dates[-1], "%Y%m%d"));
    else:
        first_epoch = [dt.datetime.strptime("19900101", "%Y%m%d") for x in name];  # placeholders
        last_epoch = [dt.datetime.strptime("20300101", "%Y%m%d") for x in name];  # placeholders
    myVelfield = Velfield(name=name, nlat=lat, elon=lon, n=n, e=e, u=u, sn=sn, se=se, su=su,
                          first_epoch=first_epoch, last_epoch=last_epoch);
    return [myVelfield];


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
    [_, dE, dN, dU, Se, Sn, Su] = np.loadtxt(filename, usecols=(2, 8, 10, 12, 14, 15, 16), skiprows=1,
                                             unpack=True);

    dtarray = [];
    ifile = open(filename);
    ifile.readline();
    for line in ifile:
        station_name = line.split()[0];
        yyMMMdd = line.split()[1];  # has format 07SEP19
        mydateobject = dt.datetime.strptime(yyMMMdd, "%y%b%d");
        dtarray.append(mydateobject);
    ifile.close();
    dN = [i * 1000.0 for i in dN];
    dE = [i * 1000.0 for i in dE];
    dU = [i * 1000.0 for i in dU];
    Sn = [i * 1000.0 for i in Sn];
    Se = [i * 1000.0 for i in Se];
    Su = [i * 1000.0 for i in Su];

    [lon, lat] = get_coordinates_for_stations([station_name], coordinates_file);  # format [lat, lon]
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


def read_pbo_hydro_file(filename, coords_file=None):
    # Useful for reading hydrology files like NLDAS, GLDAS, etc.
    # In the normal pipeline for this function, it is guaranteed to be given a real file.
    print("Reading %s" % filename);
    dtarray = [];
    station_name = filename.split('/')[-1][0:4];
    station_name = station_name.upper();
    if coords_file is not None:
        [lon, lat] = get_coordinates_for_stations([station_name], coords_file);  # format [lat, lon]
        coords = [lon[0], lat[0]];
    else:
        coords = [None, None];  # can return an object without meaningful coordinates if not asked for them.
    [dts, dN, dE, dU] = np.loadtxt(filename, usecols=(0, 3, 4, 5), dtype={'names': ('dts', 'dN', 'dE', 'dU'),
                                                                          'formats': (
                                                                              'U10', np.float, np.float, np.float)},
                                   skiprows=20, delimiter=',', unpack=True);
    for i in range(len(dts)):
        dtarray.append(dt.datetime.strptime(dts[i], "%Y-%m-%d"));
    Sn = [0.2 for i in range(len(dN))];
    Se = [0.2 for i in range(len(dE))];
    Su = [0.2 for i in range(len(dU))];

    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    return [myData];


def read_lsdm_file(filename):
    # Useful for reading hydrology files from LSDM German loading product
    # In the normal pipeline for this function, it is guaranteed to be given a real file.
    print("Reading %s" % filename);
    dtarray = [];
    station_name = filename.split('/')[-1][0:4];
    [lon, lat] = get_coordinates_for_stations([station_name]);  # format [lat, lon]
    [dts, dU, dN, dE] = np.loadtxt(filename, usecols=(0, 1, 2, 3), dtype={'names': ('dts', 'dN', 'dE', 'dU'),
                                                                          'formats': (
                                                                              'U10', np.float, np.float, np.float)},
                                   skiprows=3, delimiter=',', unpack=True);
    for i in range(len(dts)):
        dtarray.append(dt.datetime.strptime(dts[i][0:10], "%Y-%m-%d"));
    dN = [i * 1000.0 for i in dN];
    dE = [i * 1000.0 for i in dE];
    dU = [i * 1000.0 for i in dU];
    Sn = [0.2 for i in range(len(dN))];
    Se = [0.2 for i in range(len(dE))];
    Su = [0.2 for i in range(len(dU))];
    coords = [lon[0], lat[0]];
    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    return [myData];


def get_coordinates_for_stations(station_names, coordinates_file):
    # station_names is an array
    lon = [];
    lat = [];
    reference_names = [];
    reference_lons = [];
    reference_lats = [];

    # Read the file
    ifile = open(coordinates_file, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "#":
            continue;
        reference_names.append(temp[0]);
        reference_lats.append(float(temp[1]));
        testlon = float(temp[2]);
        if testlon > 180:
            testlon = testlon - 360;
        reference_lons.append(testlon);
    ifile.close();

    # find the stations
    for i in range(len(station_names)):
        myindex = reference_names.index(station_names[i]);
        lon.append(reference_lons[myindex]);
        lat.append(reference_lats[myindex]);
        if myindex == []:
            print("Error! Could not find coordinates for station %s " % station_names[i]);
            print("Returning [0,0]. ");
            lon.append(0.0);
            lat.append(0.0);

    return [lon, lat];


def get_start_times_for_unr_stations(station_names, coordinates_file):
    # Meant for UNR stations
    end_time = [];
    start_time = [];
    reference_names = [];
    reference_start_time = [];
    reference_end_time = [];

    # Read the file
    ifile = open(coordinates_file, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "#":
            continue;
        reference_names.append(temp[0]);
        reference_start_time.append(temp[7]);
        reference_end_time.append(temp[8]);
    ifile.close();

    # find the stations
    for i in range(len(station_names)):
        myindex = reference_names.index(station_names[i]);
        start_time.append(dt.datetime.strptime(reference_start_time[myindex], '%Y-%m-%d'));
        end_time.append(dt.datetime.strptime(reference_end_time[myindex], '%Y-%m-%d'));
        if myindex == []:
            print("Error! Could not find startdate for station %s " % station_names[i]);
            print("Returning [0,0]. ");
            start_time.append(dt.datetime.strptime('2000-01-01', '%Y-%m-%d'));
            end_time.append(dt.datetime.strptime('2000-01-01', '%Y-%m-%d'));
    return [start_time, end_time];


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
        "Format: lon(deg) lat(deg) e(mm) n(mm) u(mm) Se(mm) Sn(mm) Su(mm) first_date(yyyymmdd) last_date(yyyymmdd) name\n");
    for i in range(len(myVelfield.name)):
        first_epoch = dt.datetime.strftime(myVelfield.first_epoch[i], '%Y%m%d');
        last_epoch = dt.datetime.strftime(myVelfield.last_epoch[i], '%Y%m%d');
        ofile.write("%f %f %f %f %f %f %f %f %s %s %s\n" % (
            myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.u[i], myVelfield.se[i],
            myVelfield.sn[i], myVelfield.su[i], first_epoch, last_epoch, myVelfield.name[i]));
    ofile.close();
    return;


def write_gmt_velfile(myVelfield, outfile):
    ofile = open(outfile, 'w');
    ofile.write("# Format: lon(deg) lat(deg) e(mm) n(mm) Se(mm) Sn(mm) 0 0 1 name\n");
    for i in range(len(myVelfield.name)):
        if myVelfield.sn[i] < 0.2:  # trying to make a clean dataset
            ofile.write("%f %f %f %f %f %f 0 0 1 %s\n" % (
                myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i],
                myVelfield.sn[i], myVelfield.name[i]));
    ofile.close();
    return;
