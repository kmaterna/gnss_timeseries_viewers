
"""
File to read and write data from USGS formats
"""
import datetime as dt
import sys

import numpy as np
from gps_tools.gps_objects import Station_Vel, Timeseries


def usgs_vel_file_from_tsfile(infile):
    """
    Network parsing and plumbing
    We assume a parallel directory above the TS directory with a bunch of Velocity files
    From which we can derive coordinates
    """
    usgs_network = infile.split('/')[-2];
    usgs_directory = '';
    for i in range(len(infile.split('/')) - 3):
        usgs_directory = usgs_directory + infile.split('/')[i];
        usgs_directory = usgs_directory + '/'
    network_vel_file = usgs_directory + 'Velocities/' + usgs_network + '/NAM_' + usgs_network + '_vels.txt';
    return network_vel_file;


def usgs_network_from_velfile(velfile):
    """
    velfile is something like 'ITRF_Pacific_Northwest_vels.txt'
    """
    usgs_network = velfile.split('/')[-1][0:-9];
    if usgs_network[0:4] == 'NAM_':
        usgs_network = usgs_network[4:];
        usgs_refframe = 'NAM';
    elif usgs_network[0:5] == 'ITRF_':
        usgs_network = usgs_network[5:];
        usgs_refframe = 'ITRF';
    else:
        usgs_network, usgs_refframe = '', '';
    survey_flag = '_SGPS' in velfile;
    return usgs_network, usgs_refframe, survey_flag;


def read_usgs_velfile(infile, cache_file):
    """
    Read a USGS velocity file. Requires a cache of start and end dates.
    """
    print("Reading %s" % infile);
    usgs_sub_network, usgs_refframe, survey_flag = usgs_network_from_velfile(infile);
    [names, lon, lat, e, n, se, sn, u, su] = np.loadtxt(infile, skiprows=3, usecols=(0, 1, 2, 4, 5, 6, 7, 9, 10),
                                                        unpack=True, dtype={
            'names': ('name', 'lon', 'lat', 'evel', 'nvel', 'se', 'sn', 'u', 'su'),
            'formats': ('U4', float, float, float, float, float, float, float, float)})
    # Populating the first_epoch and last_epoch with information from the associated cache.
    [cache_names, startdate, enddate, subnetwork] = np.loadtxt(cache_file, unpack=True, usecols=(0, 3, 4, 5),
                                                               dtype={
                                                                   'names': ('name', 'startdate', 'enddate', 'subnet'),
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
                    print("ERROR! Station %s has no first_epoch! Stopping." % station_name);
                    sys.exit(0);
                if last_epoch is None:
                    print("ERROR! Station %s has no last_epoch! Stopping." % station_name);
                    sys.exit(0);

        myStationVel = Station_Vel(name=names[i], elon=lon[i], nlat=lat[i], e=e[i], n=n[i], u=u[i],
                                   se=se[i], sn=sn[i], su=su[i], first_epoch=first_epoch, last_epoch=last_epoch,
                                   proccenter='usgs', subnetwork=usgs_sub_network, refframe=usgs_refframe,
                                   survey=survey_flag, meas_type='gnss');
        myVelField.append(myStationVel);
    return [myVelField];


def read_USGS_ts_file(filename):
    print("Reading %s" % filename);
    station_name = filename.split('/')[-1][0:4].upper();  # the first four characters of the filename
    [datestrs, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, unpack=True, usecols=(0, 2, 3, 4, 6, 7, 8),
                                                    dtype={'names': ('datestrs', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su'),
                                                           'formats': ('U8',
                                                                       float, float, float, float, float, float)})
    dtarray = [dt.datetime.strptime(x, "%Y%m%d") for x in datestrs];
    vel_file = usgs_vel_file_from_tsfile(filename);  # get the associated velocity file for that USGS sub-network
    [names, lon, lat] = np.loadtxt(vel_file, skiprows=3, usecols=(0, 1, 2),
                                   unpack=True,
                                   dtype={'names': ('name', 'lon', 'lat'), 'formats': ('U4', float, float)})
    coords = None;  # default value in case station wasn't found
    for i in range(len(names)):
        if names[i] == station_name:
            coords = [lon[i], lat[i]];
    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su,
                        EQtimes=[]);
    return [myData];