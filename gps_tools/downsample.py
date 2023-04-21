"""
Tools to down-sample GNSS time series in time
"""
import gps_tools.gps_ts_functions
from . import gps_ts_functions
import matplotlib.pyplot as plt
import datetime as dt


def subsample_ts_at_points(Data0, dtlist, window_days=30, Se_default=None, Sn_default=None, Su_default=None):
    """
    Return a ts-object with a smaller number of time series points.

    :param Data0: GPS timeseries object
    :param dtlist: list of datetimes for downsampling
    :param window_days: integer
    :param Se_default: float
    :param Sn_default: float
    :param Su_default: float
    :return: GPS timeseries object
    """
    dE, dN, dU, Se, Sn, Su = [], [], [], [], [], [];
    for day in dtlist:
        sample_pos = gps_ts_functions.subsample_in_time(Data0, day, window_days);
        dE.append(sample_pos[0]);
        dN.append(sample_pos[1]);
        dU.append(sample_pos[2]);
        Se.append(Se_default);
        Sn.append(Sn_default);
        Su.append(Su_default);
    downsampled_object = gps_tools.gps_ts_functions.Timeseries(name=Data0.name, coords=Data0.coords, dtarray=dtlist,
                                                               dE=dE, dN=dN,
                                                               dU=dU, Sn=Sn, Se=Se, Su=Su, EQtimes=Data0.EQtimes);
    return downsampled_object;


def plot_subsampled_ts(full_station, ds_station, outname, buffer=1035):
    """
    Plot the full time series and the downsampled time series above, for visual inspection.

    :param full_station: a timeseries object
    :param ds_station: a downsampled timeseries object
    :param outname: string, filename
    :param buffer: number of days on each side to plot
    :return:
    """
    startlim = ds_station.dtarray[0] - dt.timedelta(days=buffer);
    endlim = ds_station.dtarray[-1] + dt.timedelta(days=buffer);

    f, axarr = plt.subplots(3, 1, figsize=(12, 8), dpi=300);
    axarr[0].plot(full_station.dtarray, full_station.dE, '.');
    axarr[0].set_xlim([startlim, endlim]);
    axarr[0].set_ylabel('East (mm)')
    for i in range(len(ds_station.dtarray)):
        axarr[0].plot(ds_station.dtarray[i], ds_station.dE[i], '.', color='red', markersize=15);
    axarr[1].plot(full_station.dtarray, full_station.dN, '.');
    axarr[1].set_xlim([startlim, endlim]);
    axarr[1].set_ylabel('North (mm)');
    for i in range(len(ds_station.dtarray)):
        axarr[1].plot(ds_station.dtarray[i], ds_station.dN[i], '.', color='red', markersize=15);
    axarr[2].plot(full_station.dtarray, full_station.dU, '.');
    axarr[2].set_xlim([startlim, endlim]);
    axarr[2].set_ylabel('Up (mm)');
    for i in range(len(ds_station.dtarray)):
        axarr[2].plot(ds_station.dtarray[i], ds_station.dU[i], '.', color='red', markersize=15);
    if len(ds_station.dtarray) == 2:
        axarr[0].plot([ds_station.dtarray[0], ds_station.dtarray[1]], [ds_station.dE[0], ds_station.dE[1]], color='red')
        axarr[1].plot([ds_station.dtarray[0], ds_station.dtarray[1]], [ds_station.dN[0], ds_station.dN[1]], color='red')
        axarr[2].plot([ds_station.dtarray[0], ds_station.dtarray[1]], [ds_station.dU[0], ds_station.dU[1]], color='red')
    plt.savefig(outname);
    plt.close();
    return;
