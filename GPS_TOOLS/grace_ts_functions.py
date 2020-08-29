# GRACE FUNCTIONS
import numpy as np
import matplotlib.pyplot as plt
import collections
import datetime as dt
import gps_ts_functions

Timeseries = collections.namedtuple("Timeseries", [
    'name', 'coords', 'dtarray',
    'dN', 'dE', 'dU',
    'Sn', 'Se', 'Su', 'EQtimes']);  # in mm
GRACE_TS_Array = collections.namedtuple('GRACE_TS', [
    'name', 'coords',
    'decyear', 'dtarray',
    'u', 'v', 'w'])
Paired_TS = collections.namedtuple('Paired_TS', [
    'dtarray',
    'north', 'east', 'vert',
    'N_err', 'E_err', 'V_err',
    'u', 'v', 'w']);


def input_GRACE_individual_station(filename):
    # THE GRACE DATA
    station_name = filename.split('/')[-1];  # this is the local name of the file
    station_name = station_name.split('_')[1];  # this is the four-character name
    try:
        [lon, lat, temp, u, v, w] = np.loadtxt(filename, usecols=range(1, 7), unpack=True);
    except FileNotFoundError:
        print("ERROR! Cannot find GRACE model for file %s" % filename);

    u = np.array(u);
    v = np.array(v);
    w = np.array(w);
    grace_t = get_grace_datetimes(filename);
    grace_t = [i + dt.timedelta(days=15) for i in
               grace_t];  # we add 15 days to plot the GRACE data at the center of the bin.
    grace_decyear = gps_ts_functions.get_float_times(grace_t);  # the decimal years of all the grace obs points.
    myGRACE_TS = GRACE_TS_Array(name=station_name, coords=[lon[0], lat[0]], decyear=grace_decyear, dtarray=grace_t, u=u,
                                v=v, w=w);
    return myGRACE_TS;


def pair_GPSGRACE(GPS_TS, GRACE_TS):
    # This resamples the GRACE data to match GPS that is within the range of GRACE, and forms a common time axis.
    gps_decyear = gps_ts_functions.get_float_times(GPS_TS.dtarray)
    decyear = [];
    dtarray = [];
    north_gps = [];
    east_gps = [];
    vert_gps = [];
    N_err = [];
    E_err = [];
    V_err = [];
    for i in range(len(GPS_TS.dtarray)):  # this if-statement is happening because GPS is more current than GRACE
        if min(GRACE_TS.dtarray) < GPS_TS.dtarray[i] < max(GRACE_TS.dtarray):
            decyear.append(gps_decyear[i]);
            dtarray.append(GPS_TS.dtarray[i])
            north_gps.append(GPS_TS.dN[i]);
            east_gps.append(GPS_TS.dE[i]);
            vert_gps.append(GPS_TS.dU[i]);
            N_err.append(GPS_TS.Sn[i]);
            E_err.append(GPS_TS.Se[i]);
            V_err.append(GPS_TS.Su[i]);
    grace_u = np.interp(decyear, GRACE_TS.decyear, GRACE_TS.u);
    grace_v = np.interp(decyear, GRACE_TS.decyear, GRACE_TS.v);
    grace_w = np.interp(decyear, GRACE_TS.decyear, GRACE_TS.w);
    my_paired_ts = Paired_TS(dtarray=dtarray, north=north_gps, east=east_gps, vert=vert_gps, N_err=N_err, E_err=E_err,
                             V_err=V_err, u=grace_u, v=grace_v, w=grace_w);
    return my_paired_ts;


def get_grace_datetimes(tsfile):
    ifile = open(tsfile, 'r');
    dateobjects = [];
    for line in ifile:
        temp = line.split();
        raw_string = temp[0];  # This is in the format 01-Jan-2012_31-Jan-2012
        datestring = raw_string[0:11];
        myobject = dt.datetime.strptime(datestring, '%d-%b-%Y');
        dateobjects.append(myobject);
    ifile.close();
    return dateobjects;


def get_slope(Data0, starttime=None, endtime=None):
    # Model the data with a best-fit y = mx + b.
    if starttime is None:
        starttime = Data0.dtarray[0];
    if endtime is None:
        endtime = Data0.dtarray[-1];

    # Defensive programming
    if starttime < Data0.dtarray[0]:
        starttime = Data0.dtarray[0];
    if endtime > Data0.dtarray[-1]:
        endttime = Data0.dtarray[-1];
    if endtime < Data0.dtarray[0]:
        print("Error: end time before start of array for station %s. Returning Nan" % Data0.name);
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan];
    if starttime > Data0.dtarray[-1]:
        print("Error: start time after end of array for station %s. Returning Nan" % Data0.name);
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan];

    # Cut to desired window, and remove nans
    mydtarray = [];
    myeast = [];
    mynorth = [];
    myup = [];
    for i in range(len(Data0.dtarray)):
        if starttime <= Data0.dtarray[i] <= endtime and ~np.isnan(Data0.u[i]):
            mydtarray.append(Data0.dtarray[i]);
            myeast.append(Data0.u[i]);
            mynorth.append(Data0.v[i]);
            myup.append(Data0.w[i]);

    if len(mydtarray) == 0:
        print("Error: no time array for station %s. Returning Nan" % Data0.name);
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan];
    time_duration = mydtarray[-1] - mydtarray[0];
    if time_duration.days < 365:
        print(
            "Error: using less than one year of data to estimate parameters for station %s. Returning Nan" % Data0.name);
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan];

    # doing the inversion here, since it's only one line.
    decyear = gps_ts_functions.get_float_times(mydtarray);
    east_coef = np.polyfit(decyear, myeast, 1);
    north_coef = np.polyfit(decyear, mynorth, 1);
    vert_coef = np.polyfit(decyear, myup, 1);
    east_slope = east_coef[0];
    north_slope = north_coef[0];
    vert_slope = vert_coef[0];

    # How bad is the fit to the line?
    east_trend = [east_coef[0] * x + east_coef[1] for x in decyear];
    east_detrended = [myeast[i] - east_trend[i] for i in range(len(myeast))];
    east_std = np.std(east_detrended);
    north_trend = [north_coef[0] * x + north_coef[1] for x in decyear];
    north_detrended = [mynorth[i] - north_trend[i] for i in range(len(mynorth))];
    north_std = np.std(north_detrended);
    vert_trend = [vert_coef[0] * x + vert_coef[1] for x in decyear];
    vert_detrended = [myup[i] - vert_trend[i] for i in range(len(myup))];
    vert_std = np.std(vert_detrended);

    return [east_slope, north_slope, vert_slope, east_std, north_std, vert_std];


def get_linear_annual_semiannual(Data0, starttime=None, endtime=None):
    # Model the data with a best-fit GPS = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt) + E*t + F;
    if starttime is None:
        starttime = Data0.dtarray[0];
    if endtime is None:
        endtime = Data0.dtarray[-1];

    # Defensive programming
    if starttime < Data0.dtarray[0]:
        starttime = Data0.dtarray[0];
    if endtime > Data0.dtarray[-1]:
        endtime = Data0.dtarray[-1];
    if endtime < Data0.dtarray[0]:
        print("Error: end time before start of array for station %s. Returning Nan" % Data0.name);
        east_params = [np.nan, 0, 0, 0, 0];
        north_params = [np.nan, 0, 0, 0, 0];
        up_params = [np.nan, 0, 0, 0, 0];
    if starttime > Data0.dtarray[-1]:
        print("Error: start time after end of array for station %s. Returning Nan" % Data0.name);
        east_params = [np.nan, 0, 0, 0, 0];
        north_params = [np.nan, 0, 0, 0, 0];
        up_params = [np.nan, 0, 0, 0, 0];

    # Cut to desired time window, and remove nans.
    mydtarray = [];
    myeast = [];
    mynorth = [];
    myup = [];
    for i in range(len(Data0.dtarray)):
        if starttime <= Data0.dtarray[i] <= endtime and ~np.isnan(Data0.u[i]):
            mydtarray.append(Data0.dtarray[i]);
            myeast.append(Data0.u[i]);
            mynorth.append(Data0.v[i]);
            myup.append(Data0.w[i]);

    if len(mydtarray) < 365 / 30:
        print(
            "Error: using less than one year of data to estimate parameters for station %s. Returning Nan" % Data0.name);
        east_params = [np.nan, 0, 0, 0, 0];
        north_params = [np.nan, 0, 0, 0, 0];
        up_params = [np.nan, 0, 0, 0, 0];
        return [east_params, north_params, up_params];

    decyear = gps_ts_functions.get_float_times(mydtarray);
    east_params_unordered = gps_ts_functions.invert_linear_annual_semiannual(decyear, myeast);
    north_params_unordered = gps_ts_functions.invert_linear_annual_semiannual(decyear, mynorth);
    vert_params_unordered = gps_ts_functions.invert_linear_annual_semiannual(decyear, myup);

    # The definition for returning parameters:
    # slope, a2(cos), a1(sin), s2, s1.
    east_params = [east_params_unordered[4], east_params_unordered[0], east_params_unordered[1],
                   east_params_unordered[2], east_params_unordered[3]];
    north_params = [north_params_unordered[4], north_params_unordered[0], north_params_unordered[1],
                    north_params_unordered[2], north_params_unordered[3]];
    vert_params = [vert_params_unordered[4], vert_params_unordered[0], vert_params_unordered[1],
                   vert_params_unordered[2], vert_params_unordered[3]];

    return [east_params, north_params, vert_params];


def plot_grace(station_name, filename, out_dir):
    grace_ts = input_GRACE_individual_station(filename);

    plt.figure();
    plt.plot_date(grace_ts.dtarray, grace_ts.u, '-b');
    plt.plot_date(grace_ts.dtarray, grace_ts.v, '-g');
    plt.plot_date(grace_ts.dtarray, grace_ts.w, '-r');
    plt.legend(['east', 'north', 'vertical']);
    plt.grid(True);
    plt.xlabel('Time');
    plt.ylabel('Displacement (mm)');
    plt.savefig(out_dir + station_name + "_gracets.eps");

    return;
