
"""
File to read and write data from USGS formats
"""
import datetime as dt
import glob, os, sys, subprocess
import numpy as np
from ..vel_functions import Station_Vel
from ..gps_ts_functions import Timeseries
from ..offsets import Offset


def usgs_vel_file_from_tsfile(infile):
    """
    Network parsing and plumbing
    We assume a parallel directory above TS directory with a bunch of Velocity files, from which we derive coordinates
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

    :param velfile: string, filename
    :returns: string, string, bool
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

    :param infile: string, filename
    :param cache_file: string, filename
    :returns: list of StationVel objects
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
    return myVelField;


def read_USGS_ts_file(filename) -> Timeseries:
    """
    :param filename: string, name of .rneu file for GNSS time series
    :returns: TimeSeries object
    """
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
    return myData;


def parse_offset_table_usgs(table, offset_type):
    """
    Separate a table (from a grep result) into lists of antenna and maintenance Offsets for USGS data

    :param table: string, block of text
    :param offset_type: string, either 'earthquake' or 'antenna/other'
    :return: list of Offset objects
    """
    offset_list = [];
    if len(table) == 0:
        return [];
    tablesplit = table.split('\n');
    for item in tablesplit:  # for each earthquake or offset
        if offset_type == 'earthquake':
            if 'earthquake' in item:
                words = item.split();
                if len(words) < 8:
                    continue;
                evdt = dt.datetime.strptime(words[1], "%Y-%m-%d");
                n_offset, e_offset, u_offset = float(words[3]), float(words[5]), float(words[7])
                offi = Offset(e_offset=e_offset, n_offset=n_offset, u_offset=u_offset, evdt=evdt);
                offset_list.append(offi);
        else:
            if 'earthquake' not in item:
                words = item.split();
                if len(words) < 8:
                    continue;
                evdt = dt.datetime.strptime(words[1], "%Y-%m-%d");
                n_offset, e_offset, u_offset = float(words[3]), float(words[5]), float(words[7])
                offi = Offset(e_offset=e_offset, n_offset=n_offset, u_offset=u_offset, evdt=evdt);
                offset_list.append(offi);
    return offset_list;


def search_file_for_usgs_offsets(station_name, filename):
    """
    :param station_name: string
    :param filename: string, a filename
    :returns: string, a block of text
    """
    try:
        table = subprocess.check_output("grep " + station_name + " " + filename, shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    return table;


def query_usgs_network_name(station_name, gps_ts_dir):
    """
    Given that USGS puts stations into networks, which network is a given station a member of?
    This function may print more than one.

    :param station_name: string
    :param gps_ts_dir: string, path to directory that contains subnetworks
    :returns: list of strings
    """
    network_list = [];
    directories = glob.glob(gps_ts_dir+"/*");
    print("Querying for station %s in USGS sub-networks... " % station_name);
    for item in directories:
        if os.path.isfile(item+'/'+station_name.lower()+'_NAfixed.rneu'):
            print("Found %s in %s" % (station_name, item) );
            network_list.append(item);
    print("")
    return network_list;

def read_usgs_highrate_file(filename):
    """  time format 2023/05/10 00:05:00
    :param filename: string, name of file
    :return: Timeseries object
    """
    print("Reading %s" % filename);
    dta, N, E, U = np.loadtxt(filename, usecols=(0, 2, 3, 4), delimiter=',', unpack=True, skiprows=1,
                              dtype={'names': ('dts', 'dN', 'dE', 'dU'), 'formats': ('U20', float, float, float)});
    dtarray = [dt.datetime.strptime(x, '%Y/%m/%d %H:%M:%S') for x in dta];
    myData = Timeseries(name="", coords=None, dtarray=dtarray,
                        dN=np.multiply(N, 1000), dE=np.multiply(E, 1000), dU=np.multiply(U, 1000),
                        Sn=[], Se=[], Su=[], EQtimes=[]);
    return myData;
