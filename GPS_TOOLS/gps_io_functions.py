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
Params = collections.namedtuple("Params", ['pbo_gps_dir', 'unr_gps_dir', 'pbo_earthquakes_dir', 'pbo_offsets_dir',
                                           'unr_offsets_dir', 'unr_coords_file', 'velocity_dir', 'pbo_velocities',
                                           'unr_velocities', 'gldas_dir', 'nldas_dir', 'noah_dir', 'grace_dir',
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
    pbo_gps_dir = config.get('py-config', 'pbo_gps_dir');
    unr_gps_dir = config.get('py-config', 'unr_gps_dir');
    pbo_earthquakes_dir = config.get('py-config', 'pbo_earthquakes_dir');
    unr_offsets_dir = config.get('py-config', 'unr_offsets_dir');
    pbo_offsets_dir = config.get('py-config', 'pbo_offsets_dir');
    unr_coords_file = config.get('py-config', 'unr_coords_file');
    velocity_dir = config.get('py-config', 'velocity_dir');
    pbo_velocities = config.get('py-config', 'pbo_velocities');
    unr_velocities = config.get('py-config', 'unr_velocities');
    blacklist = config.get('py-config', 'blacklist');
    gldas_dir = config.get('py-config', 'gldas_dir');
    nldas_dir = config.get('py-config', 'nldas_dir');
    noah_dir = config.get('py-config', 'noah_dir');
    grace_dir = config.get('py-config', 'grace_dir');
    lsdm_dir = config.get('py-config', 'lsdm_dir');
    stl_dir = config.get('py-config', 'stl_dir');

    myParams = Params(pbo_gps_dir=pbo_gps_dir, unr_gps_dir=unr_gps_dir, pbo_earthquakes_dir=pbo_earthquakes_dir,
                      pbo_offsets_dir=pbo_offsets_dir, unr_offsets_dir=unr_offsets_dir, unr_coords_file=unr_coords_file,
                      velocity_dir=velocity_dir, pbo_velocities=pbo_velocities,
                      unr_velocities=unr_velocities, gldas_dir=gldas_dir, nldas_dir=nldas_dir, noah_dir=noah_dir,
                      grace_dir=grace_dir, lsdm_dir=lsdm_dir, stl_dir=stl_dir, blacklist=blacklist);
    return myParams;


def read_pbo_vel_file(infile):
    # Meant for reading velocity files from the PBO/UNAVCO website.
    # Returns a Velfield object.
    print("Reading %s" % infile);
    start = 0;
    ifile = open(infile, 'r');
    name = [];
    nlat = [];
    elon = [];
    n = [];
    e = [];
    u = [];
    sn = [];
    se = [];
    su = [];
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
    n = [];
    e = [];
    u = [];
    sn = [];
    se = [];
    su = [];
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
    [first_epoch, last_epoch] = get_start_times_for_stations(name, coordinate_file);

    myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch,
                          last_epoch=last_epoch);
    return [myVelfield];


def clean_velfield(velfield, num_years=0, max_sigma=1000, coord_box=(-180, 180, -90, 90)):
    # Take the raw GPS velocities, and clean them up.
    # Remove data that's less than num_years long,
    # has formal uncertainties above max_sigma,
    # or is outside our box of interest.
    # Default arguments are meant to be global.
    name = [];
    nlat = [];
    elon = [];
    n = [];
    e = [];
    u = [];
    sn = [];
    se = [];
    su = [];
    first_epoch = [];
    last_epoch = [];
    for i in range(len(velfield.n)):
        if velfield.sn[i] > max_sigma:  # too high sigma, please exclude
            continue;
        if velfield.se[i] > max_sigma:
            continue;
        deltatime = velfield.last_epoch[i] - velfield.first_epoch[i];
        if deltatime.days <= num_years * 365.24:  # too short time interval, please exclude
            continue;
        if coord_box[0] < velfield.elon[i] < coord_box[1] and coord_box[2] < velfield.nlat[i] < coord_box[3]:
            # The station is within the box of interest.
            name.append(velfield.name[i]);
            nlat.append(velfield.nlat[i]);
            elon.append(velfield.elon[i]);
            n.append(velfield.n[i]);
            sn.append(velfield.sn[i]);
            e.append(velfield.e[i]);
            se.append(velfield.se[i]);
            u.append(velfield.u[i]);
            su.append(velfield.su[i]);
            first_epoch.append(velfield.first_epoch[i]);
            last_epoch.append(velfield.last_epoch[i]);
    myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch,
                          last_epoch=last_epoch);
    return [myVelfield];


def remove_duplicates(velfield):
    name = [];
    nlat = [];
    elon = [];
    n = [];
    e = [];
    u = [];
    sn = [];
    se = [];
    su = [];
    first_epoch = [];
    last_epoch = [];

    for i in range(len(velfield.n)):
        is_duplicate = 0;
        for j in range(len(name)):
            if abs(nlat[j] - velfield.nlat[i]) < 0.0005 and abs(elon[j] - velfield.elon[i]) < 0.0005:
                # we found a duplicate measurement.
                is_duplicate = 1;
            # Right now assuming all entries at the same lat/lon have the same velocity values.

        if is_duplicate == 0:
            name.append(velfield.name[i]);
            nlat.append(velfield.nlat[i]);
            elon.append(velfield.elon[i]);
            n.append(velfield.n[i]);
            sn.append(velfield.sn[i]);
            e.append(velfield.e[i]);
            se.append(velfield.se[i]);
            u.append(velfield.u[i]);
            su.append(velfield.su[i]);
            first_epoch.append(velfield.first_epoch[i]);
            last_epoch.append(velfield.last_epoch[i]);

    myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch,
                          last_epoch=last_epoch);
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


def get_start_times_for_stations(station_names, coordinates_file):
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
