"""
File to read and write data from Creepmeter formats
"""

import numpy as np
import pandas
import datetime as dt
from ..gps_ts_functions import Timeseries


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


def read_creepmeter_ross_road_three_col(filename) -> Timeseries:
    """
    Read the Ross Road timeseries data in a text file

    :param filename: string, filename
    :return: a Timeseries object
    """
    print("Reading file %s " % filename)
    [datestrs, num1, num2] = np.loadtxt(filename, skiprows=1, delimiter=',', usecols=(0, 1, 2), unpack=True,
                                        dtype={'names': ('dts', 'num1', 'num2'),
                                               'formats': ('U19', float, float)})
    zeros = np.zeros(np.shape(num1))
    dtarray = []
    for item in datestrs:
        dtarray.append(dt.datetime.strptime(item, "%m/%d/%Y %H:%M:%S"))
    creepmeter_data = Timeseries(name='', coords=[], dtarray=dtarray, dE=num2, dN=zeros, dU=zeros,
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
