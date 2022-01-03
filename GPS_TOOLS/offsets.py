"""
Toolbox that operates on Timeseries objects to deal with antenna offsets and earthquake offsets.
"""

import numpy as np
import collections
import datetime as dt
from . import gps_io_functions

# The namedtuple definition.  Offsets should be in mm. One object per offset
Offsets = collections.namedtuple("Offsets", ['e_offsets', 'n_offsets', 'u_offsets', 'evdts']);


def remove_offsets(Data0, offsets_obj):
    """
    the actual subtraction of offsets.
    """
    if len(offsets_obj) == 0:
        return Data0;
    newN, newE, newU = [], [], [];

    # Removing offsets
    for i in range(len(Data0.dtarray)):
        # For each day...
        tempE, tempN, tempU = Data0.dE[i], Data0.dN[i], Data0.dU[i];
        for item in offsets_obj:
            # print("removing %f mm from east at %s" % (item.e_offsets, item.evdts));
            if Data0.dtarray[i] == item.evdts:  # removing the date of the offset directly (can be messed up)
                tempE, tempN, tempU = np.nan, np.nan, np.nan;
            if Data0.dtarray[i] > item.evdts:
                tempE = tempE - item.e_offsets;
                tempN = tempN - item.n_offsets;
                tempU = tempU - item.u_offsets;
        newE.append(tempE);
        newN.append(tempN);
        newU.append(tempU);

    newData = gps_io_functions.Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=newN, dE=newE,
                                          dU=newU, Sn=Data0.Sn, Se=Data0.Se, Su=Data0.Su, EQtimes=Data0.EQtimes);
    return newData;


def fit_single_offset(dtarray, data, interval, offset_num_days):
    """
    Loop through a 1D array and find dates that are used for offset calculation.
    This is done for one component, like east
    Offset can be calculated at a day or an interval (day is just repeated twice.)
    offset_num_days is the averaging window before and after the offset time.
    """
    before_indeces, after_indeces = [], [];

    # Find the indeces of nearby days
    for i in range(len(dtarray)):
        deltat_start = dtarray[i] - interval[0];  # the beginning of the interval
        deltat_end = dtarray[i] - interval[1];  # the end of the interval
        if -offset_num_days <= deltat_start.days <= 0:
            before_indeces.append(i);
        if 0 <= deltat_end.days <= offset_num_days:
            after_indeces.append(i);

    # Identify the value of the offset.
    if before_indeces == [] or after_indeces == [] or len(before_indeces) == 1 or len(after_indeces) == 1:
        offset = 0;
        print("Warning: no data before or after offset at %s. Returning offset=0" % (
            dt.datetime.strftime(interval[0], "%Y-%m-%d")));
    else:
        before_mean = np.nanmean([data[x] for x in before_indeces]);
        after_mean = np.nanmean([data[x] for x in after_indeces]);
        offset = after_mean - before_mean;
        if offset == np.nan:
            print("Warning: np.nan offset found. Returning 0");
            offset = 0;
    return offset;


def solve_for_offsets(ts_object, offset_times, num_days=10):
    """
    Here we solve for all the offsets at a given time, which is necessary for UNR data.
    Offset_times is a list of datetime objects with unique dates.
    """
    print("Solving empirically for offsets at ", offset_times);
    Offset_obj = [];
    for i in range(len(offset_times)):
        e_offset = fit_single_offset(ts_object.dtarray, ts_object.dE, [offset_times[i], offset_times[i]], num_days);
        n_offset = fit_single_offset(ts_object.dtarray, ts_object.dN, [offset_times[i], offset_times[i]], num_days);
        u_offset = fit_single_offset(ts_object.dtarray, ts_object.dU, [offset_times[i], offset_times[i]], num_days);
        newobj = Offsets(e_offsets=e_offset, n_offsets=n_offset, u_offsets=u_offset, evdts=offset_times[i]);
        Offset_obj.append(newobj);
    return Offset_obj;


def print_offset_object(Offset_obj):
    for item in Offset_obj:
        print("%s: %.4f mmE, %.4f mmN, %.4f mmU" % (dt.datetime.strftime(item.evdts, "%Y-%m-%d"),
                                                    item.e_offsets, item.n_offsets, item.u_offsets));
    print("");
    return;
