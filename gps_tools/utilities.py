import datetime as dt
import numpy as np
import os, sys


def check_lon_sanity(lon):
    """
    Ensure longitude is between -180 and +180

    :param lon: float
    :return: float, longitude between -180 and 180
    """
    if lon > 180:
        lon = lon - 360.0;
    if lon < -360:
        lon = lon + 360;
    return lon;


def check_if_file_exists(filename):
    if os.path.isfile(filename):  # Determine if file is found on system. Provide helpful suggestions if not.
        print("Found file %s from database " % filename);
    else:  # If the file is not found on the system:
        print("Error!  Cannot find %s in database. Exiting immediately..." % filename);
        sys.exit(1);
    return;


def write_namedtuple(namedtuple_object, outfilename):
    """Utility function to serialize any generic namedtuple object into a file."""
    print('Writing file %s ' % outfilename);
    ofile = open(outfilename, 'w');
    for name, value in zip(namedtuple_object._fields, namedtuple_object):
        ofile.write(name);
        ofile.write(": ");
        ofile.write(str(value));
        ofile.write("\n");
    ofile.close();
    return;


def float_to_dt(float_time):
    # Example: 2014.194 --> datetime object
    fractional_year = str(1 + int(365.24 * (float_time - np.floor(float_time))));  # something like 004, 204, 321, etc.
    if len(fractional_year) == 1:
        fractional_year = '00' + fractional_year;
    elif len(fractional_year) == 2:
        fractional_year = '0' + fractional_year;
    if fractional_year == '367' or fractional_year == '366':
        fractional_year = '365';
    myyear = str(int(np.floor(float_time)));  # something like 2014
    my_date = dt.datetime.strptime(myyear + fractional_year, "%Y%j");
    return my_date;


def get_float_time(datetime_item):
    temp = datetime_item.strftime("%Y %j");
    temp = temp.split();
    year = temp[0]
    last_day_of_year = dt.datetime.strptime(year + "1231", "%Y%m%d");
    num_days_this_year = float(last_day_of_year.strftime("%Y %j").split()[1]);  # This is either 365 or 366
    floats = (float(temp[0]) + float(temp[1]) / num_days_this_year);
    return floats;


def get_float_times(datetimes):
    return [get_float_time(x) for x in datetimes];


def get_relative_times(datetimes, origin_dt):
    return [get_relative_time(x, origin_dt) for x in datetimes];


def get_relative_time(datetime_item, origin_dt):
    timedelta = datetime_item - origin_dt;
    relative_time = timedelta.days;
    return relative_time;


def reltime_to_dt(relative_time, origin_dt):
    timedelta = dt.timedelta(days=relative_time);
    my_date = origin_dt + timedelta;
    return my_date;


def yrnum2datetime(yearnums, starttime):
    """
    Take a set of dates, in decimal years past a certain date, and convert it into normal datetime objects.
    The input and output vectors will be exactly the same length
    """
    dtarray = [];
    for i in range(len(yearnums)):
        myyr = yearnums[i];
        dtarray.append(starttime + dt.timedelta(days=myyr * 365.24));
    return dtarray;


def get_daily_dtarray(starttime, endtime):
    """
    Produce an array of datetimes, spanning starttime to endtime at daily intervals

    :param starttime: datetime object
    :param endtime: datetime object
    :return: list of datetime objects
    """

    i = 0;
    dtarray = [starttime];
    while dtarray[-1] < endtime:
        i = i + 1;
        dtarray.append(starttime + dt.timedelta(days=i));
    return dtarray;
