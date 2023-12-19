
"""
File to read and write data from miscellaneous formats
"""
import datetime as dt
import numpy as np
import pandas, os
from .io_magnet_unr import get_coordinates_for_unr_stations
from ..vel_functions import Station_Vel
from ..gps_ts_functions import Timeseries


def read_lsdm_file(filename, coords_file=None) -> Timeseries:
    """
    Read hydrology files from LSDM German loading product
    Optional: to look up the coordinates of this station, then provide the name of a UNR-style metadata file.

    :param filename: string, name of file
    :param coords_file: optional, string, file where you can look up the station's coordinates
    :returns: TimeSeries object
    """
    print("Reading %s" % filename)
    dtarray = []
    station_name = os.path.split(filename)[1][0:4]
    if coords_file is not None:
        [lon, lat] = get_coordinates_for_unr_stations([station_name], coords_file)  # format [lat, lon]
        coords = [lon[0], lat[0]]
    else:
        coords = [None, None]  # can return an object without meaningful coordinates if not asked for them.
    [dts, dU, dN, dE] = np.loadtxt(filename, usecols=(0, 1, 2, 3), dtype={'names': ('dts', 'dN', 'dE', 'dU'),
                                                                          'formats': (
                                                                              'U10', float, float, float)},
                                   skiprows=3, delimiter=',', unpack=True)
    for i in range(len(dts)):
        dtarray.append(dt.datetime.strptime(dts[i][0:10], "%Y-%m-%d"))
    dN = [i * 1000.0 for i in dN]
    dE = [i * 1000.0 for i in dE]
    dU = [i * 1000.0 for i in dU]
    Se = 0.2 * np.ones(np.shape(dE))
    Sn = 0.2 * np.ones(np.shape(dN))
    Su = 0.2 * np.ones(np.shape(dU))
    myData = Timeseries(name=station_name, coords=coords, dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su)
    return myData


def read_grace(filename) -> Timeseries:
    """
    Read GRACE model into a GPS-style time series object.

    :param filename: string, name of file
    :returns: TimeSeries object
    """
    print("Reading %s" % filename)
    station_name = os.path.split(filename)[1]  # this is the local name of the file
    station_name = station_name.split('_')[1]  # this is the four-character name
    [dts, lon, lat, _, u, v, w] = np.loadtxt(filename, usecols=range(0, 7),
                                             dtype={'names': ('dts', 'lon', 'lat', 'elev', 'u', 'v', 'w'),
                                                    'formats': (
                                                        'U11', float, float, float, float, float, float)}, unpack=True)
    # In the file, datetime is format 01-Jan-2012_31-Jan-2012.
    # We take start date and add 15 days to put the GRACE at center of bin.
    # There will be only one point per month, generally speaking
    grace_t = [dt.datetime.strptime(x, "%d-%b-%Y") for x in dts]
    grace_t = [i + dt.timedelta(days=15) for i in grace_t]
    S = np.zeros(np.shape(u))
    GRACE_TS = Timeseries(name=station_name, coords=[lon[0], lat[0]], dtarray=grace_t, dE=u, dN=v, dU=w, Se=S, Sn=S,
                          Su=S)
    return GRACE_TS


def read_humanread_vel_file(infile):
    """
    Reading velfield file with format:
    lon(deg) lat(deg) e(mm) n(mm) u(mm) Se(mm) Sn(mm) Su(mm) first_dt(yyyymmdd) last_dt(yyyymmdd) name\n")

    :param infile: string, filename
    :returns: list of StationVel objects
    """
    myVelfield = []
    print("Reading %s" % infile)
    ifile = open(infile, 'r')
    for line in ifile:
        if line.split()[0] == '#':
            continue
        temp = line.split()
        new_station = Station_Vel(name=temp[10], elon=float(temp[0]), nlat=float(temp[1]), e=float(temp[2]),
                                  n=float(temp[3]), u=float(temp[4]), se=float(temp[5]), sn=float(temp[6]),
                                  su=float(temp[7]), first_epoch=dt.datetime.strptime(temp[8], "%Y%m%d"),
                                  last_epoch=dt.datetime.strptime(temp[9], "%Y%m%d"), refframe='',
                                  proccenter='', subnetwork='', survey=False, meas_type='gnss')
        myVelfield.append(new_station)
    ifile.close()
    return myVelfield


def write_humanread_vel_file(myVelfield, outfile):
    """
    :param myVelfield: list of StationVel objects
    :param outfile: string, filename
    """
    print("writing human-readable velfile %s" % outfile)
    ofile = open(outfile, 'w')
    ofile.write(
        "# Fmt: lon(deg) lat(deg) e(mm) n(mm) u(mm) Se(mm) Sn(mm) Su(mm) first_dt(yyyymmdd) last_dt(yyyymmdd) name\n")
    for station_vel in myVelfield:
        first_epoch = dt.datetime.strftime(station_vel.first_epoch, '%Y%m%d')
        last_epoch = dt.datetime.strftime(station_vel.last_epoch, '%Y%m%d')
        ofile.write("%f %f %f %f %f %f %f %f %s %s %s\n" % (
            station_vel.elon, station_vel.nlat, station_vel.e, station_vel.n, station_vel.u, station_vel.se,
            station_vel.sn, station_vel.su, first_epoch, last_epoch, station_vel.name))
    ofile.close()
    return


def write_stationvel_file(myVelfield, outfile, metadata=None):
    """
    :param myVelfield: list of StationVel objects
    :param outfile: string, filename
    :param metadata: string, optional, used for header
    """
    print("writing human-readable velfile in station-vel format, %s" % outfile)
    ofile = open(outfile, 'w')
    ofile.write(
        "# Format: lon(deg) lat(deg) VE(mm) VN(mm) VU(mm) SE(mm) SN(mm) SU(mm) name(optional)\n")
    if metadata:
        ofile.write("# derived from metadata: %s\n" % metadata)
    for station_vel in myVelfield:
        ofile.write("%f %f %f %f %f %f %f %f %s\n" % (
            station_vel.elon, station_vel.nlat, station_vel.e, station_vel.n, station_vel.u, station_vel.se,
            station_vel.sn, station_vel.su, station_vel.name))
    ofile.close()
    return


def write_gmt_velfile(myVelfield, outfile):
    """
    :param myVelfield: list of StationVel objects
    :param outfile: string, filename
    """
    print("Writing gmt velfile %s" % outfile)
    ofile = open(outfile, 'w')
    ofile.write("# Format: lon(deg) lat(deg) e(mm) n(mm) Se(mm) Sn(mm) 0 0 1 name\n")
    for station_vel in myVelfield:
        ofile.write("%f %f %f %f %f %f 0 0 1 %s\n" % (
            station_vel.elon, station_vel.nlat, station_vel.e, station_vel.n, station_vel.se,
            station_vel.sn, station_vel.name))
    ofile.close()
    return


def read_blacklist(blacklist_file) -> list:
    """
    :param blacklist_file: string, filename
    :return: list of strings
    """
    print("Reading blacklist file %s " % blacklist_file)
    with open(blacklist_file, "r") as f:
        blacklist = []
        blacklist_all = f.readlines()
        for line in blacklist_all:
            blacklist.append(line.split()[0])
    return blacklist


def write_stl(Data, filename):
    ofile = open(filename, 'w')
    for i in range(len(Data.dtarray)):
        timestamp = dt.datetime.strftime(Data.dtarray[i], "%Y%m%d")
        ofile.write("%s %f %f %f %f %f %f\n" % (
            timestamp, Data.dE[i], Data.dN[i], Data.dU[i], Data.Se[i], Data.Sn[i], Data.Su[i]))
    ofile.close()
    return


def read_lake_loading_ts(infile) -> Timeseries:
    """
    :param infile: string, filename
    :returns: a TimeSeries object
    """
    [dtstrings, u, v, w] = np.loadtxt(infile, usecols=(0, 4, 5, 6), unpack=True, dtype={
        'names': ('d', 'u', 'v', 'w'), 'formats': ('U10', np.float, np.float, np.float)})
    dtarray = [dt.datetime.strptime(x, "%Y-%m-%d") for x in dtstrings]
    S = np.zeros(np.shape(u))
    loading_defo = Timeseries(name='', coords=[], dtarray=dtarray, dE=u, dN=v, dU=w, Sn=S, Se=S, Su=S)
    return loading_defo


def read_creepmeter(filename) -> Timeseries:
    """
    Read an excel file for creep-meter timeseries such as those from the Imperial Valley.

    :param filename: string, filename
    :return: a TimeSeries object
    """
    print("Reading file %s " % filename)
    df = pandas.read_excel(filename, sheet_name=0)
    ts = df['UTC'].values
    dtarray = [pandas.to_datetime(x) for x in ts]
    displacement = df['dextral'].values.tolist()
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=displacement,
                                 dN=np.zeros(np.shape(displacement)), dU=np.zeros(np.shape(displacement)),
                                 Sn=np.zeros(np.shape(displacement)), Se=np.zeros(np.shape(displacement)),
                                 Su=np.zeros(np.shape(displacement)))
    return creepmeter_data

def read_creepmeter2(filename) -> Timeseries:
    """
    Read an excel file for creep-meter timeseries such as those from the Imperial Valley.

    :param filename: string, filename
    :return: a TimeSeries object
    """
    print("Reading file %s " % filename)
    df = pandas.read_excel(filename, sheet_name=0)
    ts = df['date'].values
    dtarray = [pandas.to_datetime(x) for x in ts]
    displacement = df['slip'].values.tolist()
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=displacement,
                                 dN=np.zeros(np.shape(displacement)), dU=np.zeros(np.shape(displacement)),
                                 Sn=(), Se=(), Su=())
    return creepmeter_data


def read_mit_hr_fig(filename_dt, filename_e, filename_n, is_six_hour=False) -> Timeseries:
    """
    :param filename_dt: string, filename of high-rate .fig file produced at MIT
    :param filename_e: string, filename of high-rate .fig file produced at MIT
    :param filename_n: string, filename of high-rate .fig file produced at MIT
    :param is_six_hour: bool
    :return: a TimeSeries object
    """
    print("Reading file %s " % filename_dt)
    dtarray = []
    with open(filename_dt, 'r') as f:
        for line in f:
            if len(line.split()) > 1:
                dtarray.append(dt.datetime.strptime(line.split()[0]+' '+line.split()[1], '%d-%b-%Y %H:%M:%S'))
    if is_six_hour:
        e, se = np.loadtxt(filename_e, unpack=True, usecols=(2, 3))
    else:
        e = np.loadtxt(filename_e, unpack=True)
        se = np.zeros(np.shape(e))
    n = np.loadtxt(filename_n, unpack=True)
    hr_data = Timeseries(name='P503', coords=[], dtarray=dtarray, dE=e, dN=n, dU=np.zeros(np.shape(n)),
                         Se=se, Sn=np.zeros(np.shape(n)), Su=np.zeros(np.shape(n)))
    return hr_data


def read_uw_kinematic(filename) -> Timeseries:
    """
    Reading the data produced by UW for every 30 minutes of GNSS high-rate time series

    :param filename: string
    :return:  Timeseries
    """
    print("Reading file %s " % filename)
    dtvals, e, n, u = np.loadtxt(filename, usecols=(0, 1, 2, 3), unpack=True)
    dtarray = []
    for item in dtvals:
        day_of_month = int(np.floor(item))
        fractional_day = item - day_of_month
        seconds_after_midnight = fractional_day * (24 * 60 * 60)
        hours_after_midnight = int(np.floor(seconds_after_midnight / (60 * 60)))
        hours_string = define_hours_string(hours_after_midnight)
        if seconds_after_midnight > hours_after_midnight * 60 * 60 + 1000:
            minutes_string = '45'
        else:
            minutes_string = '15'
        newdatetime = dt.datetime.strptime('2023-05-'+str(day_of_month) + 'T' +
                                           hours_string+':'+minutes_string+':00', '%Y-%m-%dT%H:%M:%S')
        dtarray.append(newdatetime)
    hr_data = Timeseries(name='P503', coords=[], dtarray=dtarray, dE=e, dN=n, dU=u, Se=np.zeros(np.shape(e)),
                         Sn=np.zeros(np.shape(n)), Su=np.zeros(np.shape(u)))
    return hr_data


def define_hours_string(hours_int):
    if hours_int <= 9:
        hours_string = '0' + str(hours_int)
    else:
        hours_string = str(hours_int)
    return hours_string;


def write_gmt_ts_file(ts_data, outfile):
    """Write a time series object into a format that GMT can subsequently plot."""
    print("Writing file %s " % outfile)
    with open(outfile, 'w') as ofile:
        for i in range(len(ts_data.dtarray)):
            ofile.write("%s %f %f %f\n" % (dt.datetime.strftime(ts_data.dtarray[i], "%Y-%m-%dT%H-%M-%S"),
                                           ts_data.dE[i],
                                           ts_data.dN[i],
                                           ts_data.dU[i]) )
    return
