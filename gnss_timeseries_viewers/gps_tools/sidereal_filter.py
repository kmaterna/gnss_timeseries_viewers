"""
Implement a sidereal filter for a high-rate time series.
For the GPS, that is 24 hours minus 236 seconds.
"""

import numpy as np
from .gps_ts_functions import Timeseries


def sidereal_filter(Data0: Timeseries, starttime, endtime) -> Timeseries:
    """

    :param Data0: a Timeseries object
    :param starttime: a datetime object
    :param endtime: a datetime object
    :return: a Timeseries object
    """

    filter_created = Data0.impose_time_limits(starttime, endtime)

    sidereal_interval_dt = endtime - starttime
    sampling_interval = filter_created.dtarray[1] - filter_created.dtarray[0]
    sampling_int_sec = sampling_interval.seconds
    print(f'Length of sidereal day: {sidereal_interval_dt}')
    print(f'Length of sidereal day in seconds: {sidereal_interval_dt.seconds}')
    print(f'Sampling interval in seconds: {sampling_interval.seconds}')
    beginning_of_filter = filter_created.dtarray[0]

    dE, dN, dU, Se, Sn, Su = [], [], [], [], [], []

    for i in range(len(Data0.dtarray)):
        normalized_time = Data0.dtarray[i] - beginning_of_filter
        target_filter_sample = np.floor(np.mod((normalized_time.seconds + sampling_int_sec/2),
                                               sidereal_interval_dt.seconds) / sampling_interval.seconds)
        target_filter_sample = int(target_filter_sample)
        dE.append(Data0.dE[i] - filter_created.dE[target_filter_sample])
        dN.append(Data0.dN[i] - filter_created.dN[target_filter_sample])
        dU.append(Data0.dU[i] - filter_created.dU[target_filter_sample])
        Se.append(Data0.Se[i])
        Sn.append(Data0.Sn[i])
        Su.append(Data0.Su[i])

    newData = Timeseries(name=Data0.name, coords=Data0.coords, dtarray=Data0.dtarray, dN=np.array(dN),
                         dE=np.array(dE), dU=np.array(dU), Sn=np.array(Sn), Se=np.array(Se),
                         Su=np.array(Su), EQtimes=Data0.EQtimes, Offsettimes=Data0.Offsettimes)

    return newData
