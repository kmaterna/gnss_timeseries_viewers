
"""
File to read and write data from Network of the Americas / Plate Boundary Observatory formats
"""
import datetime as dt

import numpy as np
from gps_tools import utilities
from gps_tools.file_io.io_magnet_unr import get_coordinates_for_unr_stations
from gps_tools.gps_objects import Station_Vel, Timeseries


def read_pbo_vel_file(infile):
    """
    Read velocity files from the PBO/UNAVCO website.
    Returns a list of station_vel objects.
    """
    myVelfield = [];
    print("Reading %s" % infile);
    start = 0;
    ifile = open(infile, 'r');
    for line in ifile:
        if start == 1:
            temp = line.split();
            name = temp[0];
            nlat, elon = float(temp[7]), float(temp[8]);
            elon = utilities.check_lon_sanity(elon);
            n = float(temp[19]) * 1000.0;
            e = float(temp[20]) * 1000.0;
            u = float(temp[21]) * 1000.0;
            sn = float(temp[22]) * 1000.0;
            se = float(temp[23]) * 1000.0;
            su = float(temp[24]) * 1000.0;
            t1, t2 = temp[-2], temp[-1];
            first_epoch = dt.datetime.strptime(t1[0:8], '%Y%m%d');
            last_epoch = dt.datetime.strptime(t2[0:8], '%Y%m%d');
            myStationVel = Station_Vel(name=name, elon=elon, nlat=nlat, e=e, n=n, u=u,
                                       se=se, sn=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch,
                                       proccenter='pbo', subnetwork='', refframe='ITRF', survey=0, meas_type='gnss');
            myVelfield.append(myStationVel);
        if "*" in line:
            start = 1;

    ifile.close();
    return [myVelfield];


def read_pbo_vel_file_format(infile):
    """
    Read velocity files from the PBO/UNAVCO website.
    Returns a list of station_vel objects.
    Splits the lines by character, like Fortran style
    """
    myVelfield = [];
    print("Reading %s" % infile);
    start = 0;
    ifile = open(infile, 'r');
    for line in ifile:
        if start == 1:
            name = line[1:5];
            nlat = float(line[97:111]);
            elon = float(line[112:127]);
            elon = utilities.check_lon_sanity(elon);

            n = float(line[214:223]) * 1000.0;
            e = float(line[223:231]) * 1000.0;
            u = float(line[232:241]) * 1000.0;
            sn = float(line[241:249]) * 1000.0;
            se = float(line[249:257]) * 1000.0;
            su = float(line[257:265]) * 1000.0;
            t1 = line.split()[-2];
            t2 = line.split()[-1];
            first_epoch = dt.datetime.strptime(t1[0:8], '%Y%m%d');
            last_epoch = dt.datetime.strptime(t2[0:8], '%Y%m%d');
            myStationVel = Station_Vel(name=name, elon=elon, nlat=nlat, e=e, n=n, u=u,
                                       se=se, sn=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch,
                                       proccenter='pbo', subnetwork='', refframe='ITRF', survey=0, meas_type='gnss');
            myVelfield.append(myStationVel);
        if "*" in line:
            start = 1;

    ifile.close();
    return [myVelfield];


def read_gamit_velfile(infile):
    """
    Read a velocity file for example from GAMIT processing
    A simpler format than the PBO (only 13 fields), without starttime and stoptime information.
    """
    myVelfield = [];
    print("Reading %s" % infile);
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "#" or temp[0][0] == "#":
            continue;
        else:
            elon, nlat = float(temp[0]), float(temp[1]);
            elon = utilities.check_lon_sanity(elon);
            e, n, se, sn = float(temp[2]), float(temp[3]), float(temp[6]), float(temp[7]);
            u, su = float(temp[9]), float(temp[11]);
            name = temp[12][0:4];
            first_epoch = dt.datetime.strptime("19900101", "%Y%m%d");  # placeholders
            last_epoch = dt.datetime.strptime("20300101", "%Y%m%d");  # placeholders

            myStationVel = Station_Vel(name=name, elon=elon, nlat=nlat, e=e, n=n, u=u,
                                       se=se, sn=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch,
                                       proccenter='gamit', subnetwork='', refframe='ITRF', survey=0, meas_type='gnss');
            myVelfield.append(myStationVel);
    ifile.close();
    return [myVelfield];


def read_pbo_pos_file(filename):
    print("Reading %s" % filename);
    [yyyymmdd, Nlat, Elong, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, skiprows=37, unpack=True,
                                                                 usecols=(0, 12, 13, 15, 16, 17, 18, 19, 20));
    Elong = utilities.check_lon_sanity(Elong[0]);
    dN = [i * 1000.0 for i in dN];
    dE = [i * 1000.0 for i in dE];
    dU = [i * 1000.0 for i in dU];
    Sn = [i * 1000.0 for i in Sn];
    Se = [i * 1000.0 for i in Se];
    Su = [i * 1000.0 for i in Su];
    specific_file = filename.split('/')[-1];
    dtarray = [dt.datetime.strptime(str(int(i)), "%Y%m%d") for i in yyyymmdd];
    myData = Timeseries(name=specific_file[0:4], coords=[Elong, Nlat[0]], dtarray=dtarray, dN=dN, dE=dE, dU=dU,
                        Sn=Sn, Se=Se, Su=Su, EQtimes=[]);
    print("Reading data for station %s in time range %s:%s" % (
        myData.name, dt.datetime.strftime(myData.dtarray[0], "%Y-%m-%d"),
        dt.datetime.strftime(myData.dtarray[-1], "%Y-%m-%d")));
    return [myData];


def read_pbo_hydro_file(filename, coords_file=None):
    """
    Read hydrology files like NLDAS, GLDAS, etc.
    In the normal pipeline for this function, it is guaranteed to be given a real file.
    """
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
                                                                          'formats': ('U10', float, float, float)},
                                   skiprows=20, delimiter=',', unpack=True);
    for i in range(len(dts)):
        dtarray.append(dt.datetime.strptime(dts[i], "%Y-%m-%d"));
    Se = 0.2 * np.ones(np.shape(dE));
    Sn = 0.2 * np.ones(np.shape(dN));
    Su = 0.2 * np.ones(np.shape(dU));

    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    return [myData];


def write_pbo_pos_file(ts_object, filename, comment=""):
    """
    Useful for writing common mode objects, etc.
    Opposite of read_pbo_pos_file(filename)
    """
    print("Writing pbo-posfile %s" % filename);
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


def restrict_pbo_vel_file(infile, outfile, coord_box):
    """
    Copying format of pbo velocities, let's make a restricted dataset.
    In the ideal case, this should probably be multiple functions (READ, FILTER, WRITE), or possibly not exist.
    """
    ifile = open(infile, 'r');
    ofile = open(outfile, 'w');
    start = 0;
    for line in ifile:
        if start == 1:
            nlat = float(line[97:111]);
            elon = float(line[112:127]);
            elon = utilities.check_lon_sanity(elon);
            if coord_box[0] <= elon <= coord_box[1] and coord_box[2] <= nlat <= coord_box[3]:
                ofile.write(line);
        if start == 0:
            ofile.write(line);
        if line.split()[0][0] == '*':
            start = 1;
            continue;
    ifile.close();
    ofile.close();
    return;