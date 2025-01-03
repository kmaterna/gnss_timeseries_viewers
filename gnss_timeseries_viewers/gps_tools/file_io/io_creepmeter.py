"""
File to read and write data from Creepmeter formats
"""

import numpy as np
import pandas
import datetime as dt
from ..gps_ts_functions import Timeseries


def read_creepmeter(filename) -> Timeseries:
    """
    Read an Excel file for creep-meter timeseries such as those from the Imperial Valley.

    :param filename: string, filename
    :return: a TimeSeries object
    """
    print("Reading file %s " % filename)
    df = pandas.read_excel(filename, sheet_name=0)
    ts = df['UTC'].values
    dtarray = [pandas.to_datetime(x) for x in ts]
    displacement = df['dextral'].values.tolist()
    zeros = np.zeros(np.shape(displacement))
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=displacement, dN=zeros, dU=zeros,
                                 Sn=zeros, Se=zeros, Su=zeros)
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
    zeros = np.zeros(np.shape(displacement))
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=displacement, dN=zeros, dU=zeros,
                                 Sn=zeros, Se=zeros, Su=zeros)
    return creepmeter_data


def read_creepmeter_ross_road_three_col(filename, column_num=2) -> Timeseries:
    """
    Read the Ross Road timeseries data in a text file

    :param filename: string, filename
    :param column_num: zero-index column number with data
    :return: a Timeseries object
    """
    print("Reading file %s " % filename)
    [datestrs, num1, num2] = np.loadtxt(filename, skiprows=1, delimiter=',', usecols=(0, 1, 2), unpack=True,
                                        dtype={'names': ('dts', 'num1', 'num2'),
                                               'formats': ('U19', float, float)})
    zeros = np.zeros(np.shape(num1))
    dtarray = []
    if column_num == 2:
        data = num2
    else:
        data = num1
    for item in datestrs:
        dtarray.append(dt.datetime.strptime(item, "%m/%d/%Y %H:%M:%S"))
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=data, dN=zeros, dU=zeros,
                                 Sn=zeros, Se=zeros, Su=zeros)
    return creepmeter_data


def read_creepmeter_ross_road_two_col(filename) -> Timeseries:
    """
    Read the Ross Road timeseries data in a text file

    :param filename: string, filename
    :return: a Timeseries object
    """
    print("Reading file %s " % filename)
    [datestrs, num1] = np.loadtxt(filename, skiprows=1, delimiter=',', usecols=(0, 1), unpack=True,
                                  dtype={'names': ('dts', 'num1'),
                                         'formats': ('U19', float)})
    zeros = np.zeros(np.shape(num1))
    dtarray = []
    for item in datestrs:
        dtarray.append(dt.datetime.strptime(item, "%m/%d/%Y %H:%M:%S"))
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=num1, dN=zeros, dU=zeros,
                                 Sn=zeros, Se=zeros, Su=zeros)
    return creepmeter_data


def read_creepmeter_nyland_ranch(filename) -> Timeseries:
    """
    Read a creepmeter time series from USGS, for example Nyland ranch, downloaded from here:
    https://earthquake.usgs.gov/monitoring/deformation/data/download.php

    NOTE: The conversion from float-day to datetime probably rounds incorrectly at the one-minute level.

    :param filename: string, filename
    :return: a Timeseries object
    """
    print("Reading file %s " % filename)
    [year, decimal_date, value] = np.loadtxt(filename, unpack=True, usecols=(0, 1, 2))
    zeros = np.zeros(np.shape(year))
    dtarray = []
    for i in range(len(year)):
        year_str = str(int(year[i]))
        doy = int(np.floor(decimal_date[i]))
        fractional_day = decimal_date[i] - doy   # something like 0.347000 (part of a day)
        number_of_minutes = fractional_day * (60*24)  # number of minutes past midnight
        hour = int(np.floor(number_of_minutes / 60))
        minutes = int(np.floor(number_of_minutes)) - hour*60  # minutes past the hour
        if hour < 10:
            hour_string = "0" + str(hour)
        else:
            hour_string = str(hour)
        if minutes < 10:
            minute_string = "0" + str(minutes)
        else:
            minute_string = str(minutes)
        if doy < 10:
            doy_str = "00"+str(doy)
        elif doy < 100:
            doy_str = "0"+str(doy)
        else:
            doy_str = str(doy)
        formatted_datetime = dt.datetime.strptime(str(year_str)+"-"+str(doy_str)+"-"+str(hour_string)
                                                  + "-" + str(minute_string), "%Y-%j-%H-%M")
        dtarray.append(formatted_datetime)
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=np.array(value),
                                 dN=zeros, dU=zeros, Se=zeros, Sn=zeros, Su=zeros)
    return creepmeter_data
