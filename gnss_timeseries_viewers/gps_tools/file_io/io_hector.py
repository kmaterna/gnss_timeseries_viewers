"""
Write a time series into hector format.
Assumes you have hectorp installed in your python environment.
"""

import subprocess
import datetime as dt


def datetime_to_momtime(today: dt.datetime):
    """
    Assumes your computer has Hectorp installed (can be installed on pip)

    :param today: datetime
    :return: float of modified-julian-day
    """
    year, month = dt.datetime.strftime(today, "%Y"), dt.datetime.strftime(today, "%m")
    day = dt.datetime.strftime(today, "%d")
    cmd = "date2mjd " + year + " " + month + " " + day + " 00 00 00"
    output = subprocess.check_output(cmd, shell=True, encoding="utf8")
    mjd = float(output.split("\n")[-2].split(":")[-1])
    return mjd


def write_hector_timeseries(Data, filename, data_array, offset_times=(), log_times=()):
    """
    Write the format of time series used for Hector MOM files (modified-julian-date, observation, model).
    The measurement output will be in millimeters.

    :param Data: Timeseries object.
    :param filename: string, filename of output file.
    :param data_array: list of floats, of same length as Data.dtarray
    :param offset_times: optional list of times when offsets are expected.
    :param log_times: optional list of times and time-constants, when log transients will occur (postseismic)
    """
    print("Writing %d data points to .mom data file %s " % (len(Data.dtarray), filename))
    if len(data_array) != len(Data.dtarray):
        raise ValueError("Error! Length of dtarray and provided data_array do not match.")
    offset_times = set(offset_times)
    with (open(filename, 'w') as f):
        f.write("# sampling period 1.0\n")  # for daily observations
        for offset_time in offset_times:
            mjd = datetime_to_momtime(offset_time)
            f.write("# offset %.1f\n" % mjd)
        for log_time in log_times:
            if Data.dtarray[0] <= log_time[0] <= Data.dtarray[-1]:
                mjd = datetime_to_momtime(log_time[0])
                f.write('# log %.1f %.1f\n' % (mjd, log_time[1]))
        for i in range(len(Data.dtarray)):
            mjd = datetime_to_momtime(Data.dtarray[i])
            f.write("%.1f %.4f\n" % (mjd, data_array[i]))
    return
