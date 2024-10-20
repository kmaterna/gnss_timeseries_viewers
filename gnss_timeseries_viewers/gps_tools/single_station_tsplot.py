"""
Make a basic python plot of single-station position time series with corrections optional.
"""

import matplotlib.pyplot as plt
import os
from . import gps_seasonal_removals, offsets, load_gnss
import datetime as dt


def view_single_station(station_name, data_config_file, offsets_remove=1, earthquakes_remove=0, outliers_remove=0,
                        seasonals_remove=0, starttime=None, endtime=None, outliers_def=15,
                        seasonals_type='lssq', datasource='cwu', plot_detrended=True, refframe='NA', outdir=""):
    """
    :param station_name: string
    :param data_config_file: string
    :param offsets_remove: bool, default True
    :param earthquakes_remove: bool, default False
    :param outliers_remove: bool, default False
    :param seasonals_remove: bool, default False
    :param starttime: dt.datetime object, default None
    :param endtime: dt.datetime object, default None
    :param outliers_def: mm, default 15
    :param seasonals_type: string, default 'lssq'
    :param datasource: string, default cwu
    :param plot_detrended: whether to plot the de-trended version, default true
    :param refframe: string, default NA
    :param outdir: string
    """
    db_params = config_db(station_name, datasource, refframe)
    plot_params = config_view(offsets_remove, earthquakes_remove, outliers_remove, outliers_def, seasonals_remove,
                              seasonals_type, starttime, endtime, plot_detrended)
    [myData, offset_obj, eq_obj] = input_data(station_name, data_config_file, db_params)
    [updatedData, detrended] = compute(data_config_file, myData, offset_obj, eq_obj, plot_params)
    single_ts_plot(updatedData, detrended, plot_params=plot_params, db_params=db_params, outdir=outdir)


# -------------- CONFIGURE ------------ #
def config_db(station, datasource, refframe):
    database_params = {"station": station, 'datasource': datasource, 'refframe': refframe}
    print("\n---------- VIEWING GNSS TIMESERIES ------------")
    print("Database Parameters: station = %s, dattasource = %s, refframe = %s " % (station, datasource, refframe))
    return database_params


def config_view(offsets_remove, earthquakes_remove, outliers_remove, outliers_def, seasonals_remove, seasonals_type,
                starttime, endtime, plot_detrended):
    """
    outliers_def : mm away from median filter.
    offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove, plot_detrended : booleans
    seasonals type : lssq, nldas, gldas, grace, lsdm, shasta, stl, notch
    """
    view_params = locals()  # package all provided arguments into a dictionary
    print("Viewing Parameters: earthquakes_remove = %d, outliers_remove = %d, seasonals_remove = %d"
          % (earthquakes_remove, outliers_remove, seasonals_remove))
    print("----------------------------------")
    return view_params


# ----------- INPUTS ---------------- # 
def input_data(st_name, data_config_file, db_params):
    database = load_gnss.create_station_repo(data_config_file, db_params['refframe'], db_params['datasource'])
    [myData, offset_obj, eq_obj] = database.load_station(st_name)
    return [myData, offset_obj, eq_obj]


# -------------- COMPUTE ------------ # 
def compute(data_config_file, myData, offset_obj, eq_obj, plot_params):
    if plot_params["starttime"] is None:
        plot_params["starttime"] = myData.dtarray[0]
    if plot_params["endtime"] is None:
        plot_params["endtime"] = myData.dtarray[-1]
    starttime = plot_params["starttime"]
    endtime = plot_params["endtime"]
    newData = myData.impose_time_limits(starttime, endtime)
    if plot_params["offsets_remove"]:  # Remove offsets and antenna changes
        newData = offsets.remove_offsets(newData, offset_obj)
    if plot_params["outliers_remove"]:  # Remove outliers
        outliers_def = plot_params["outliers_def"]
        newData = newData.remove_outliers(outliers_def)
    if plot_params["earthquakes_remove"]:  # Remove earthquakes
        newData = offsets.remove_offsets(newData, eq_obj)
    trend_out = gps_seasonal_removals.make_detrended_ts(newData, plot_params["seasonals_remove"],
                                                        plot_params["seasonals_type"], data_config_file)
    return [newData, trend_out]


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, detrended=None, plot_params=None, db_params=None, outdir="", title=None, savename=None,
                   buffer_days=0, label_rotation=0, detrended_label='detrended'):
    """
    :param ts_obj: a TimeSeries object
    :param detrended: another TimeSeries object, optional
    :param plot_params: dictionary, can be used to set the filename and title of the resulting plot
    :param db_params: dictionary, can be used to set the filename and title of the resulting plot
    :param outdir: string, default current directory
    :param title: string, optional, can be used to override default title
    :param savename: string, optional, can be used to override default destination
    :param buffer_days: number of days to buffer the plot on both ends, default 0
    :param label_rotation: rotation parameter for text on x-axis
    :param detrended_label: string, label on the mirrored y-axis
    """

    if outdir != '':
        os.makedirs(outdir, exist_ok=True)

    if not title:
        title, _ = get_figure_name(plot_params, db_params)
    if not savename:
        _, savename = get_figure_name(plot_params, db_params)

    savename = outdir + savename

    if not plot_params["plot_detrended"]:
        detrended = None

    # The major figure
    dpival = 500
    label_fontsize = 18
    markersize = 3.0
    grid_linewidth = 0.5
    eq_linewidth = 0.5

    # noinspection PyTypeChecker
    [_, axarr] = plt.subplots(3, 1, sharex=True, figsize=(10, 7), dpi=dpival)
    axarr[0].plot(ts_obj.dtarray, ts_obj.dE, color='blue', markeredgecolor='black', marker='.', markersize=markersize,
                  linestyle='')
    axarr[0].grid(linestyle='--', linewidth=grid_linewidth)
    axarr[0].set_ylabel('east (mm)', fontsize=label_fontsize)
    bottom, top = axarr[0].get_ylim()
    for i in range(len(ts_obj.EQtimes)):
        axarr[0].plot([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k', linewidth=eq_linewidth)
    for i in range(len(ts_obj.Offsettimes)):
        axarr[0].plot([ts_obj.Offsettimes[i], ts_obj.Offsettimes[i]], [bottom, top], '--c', linewidth=eq_linewidth)
    if detrended:
        ax1 = axarr[0].twinx()
        ax1.plot(detrended.dtarray, detrended.dE, marker='D', markersize=1.0, color='red', linestyle='')
        ax1.set_ylabel(detrended_label+' (mm)', fontsize=label_fontsize - 2, color='red')
        ax1.tick_params(labelcolor='red', labelsize=label_fontsize, axis='both')
    axarr[0].tick_params(labelsize=label_fontsize)

    axarr[1].plot(ts_obj.dtarray, ts_obj.dN, color='blue', markeredgecolor='black', marker='.', markersize=markersize,
                  linestyle='')
    axarr[1].grid(linestyle='--', linewidth=grid_linewidth)
    axarr[1].set_ylabel('north (mm)', fontsize=label_fontsize)
    bottom, top = axarr[1].get_ylim()
    for i in range(len(ts_obj.EQtimes)):
        axarr[1].plot([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k', linewidth=eq_linewidth)
    for i in range(len(ts_obj.Offsettimes)):
        axarr[1].plot([ts_obj.Offsettimes[i], ts_obj.Offsettimes[i]], [bottom, top], '--c', linewidth=eq_linewidth)
    if detrended:
        ax2 = axarr[1].twinx()
        ax2.plot(detrended.dtarray, detrended.dN, marker='D', markersize=1.0, color='red', linestyle='')
        ax2.set_ylabel(detrended_label+' (mm)', fontsize=label_fontsize - 2, color='red')
        ax2.tick_params(labelcolor='red', labelsize=label_fontsize, axis='both')
    axarr[1].tick_params(labelsize=label_fontsize)

    axarr[2].plot(ts_obj.dtarray, ts_obj.dU, color='blue', markeredgecolor='black', marker='.', markersize=markersize,
                  linestyle='')
    axarr[2].grid(linestyle='--', linewidth=grid_linewidth)
    axarr[2].set_ylabel('vertical (mm)', fontsize=label_fontsize)
    bottom, top = axarr[2].get_ylim()
    for i in range(len(ts_obj.EQtimes)):
        axarr[2].plot([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k', linewidth=eq_linewidth)
    for i in range(len(ts_obj.Offsettimes)):
        axarr[2].plot([ts_obj.Offsettimes[i], ts_obj.Offsettimes[i]], [bottom, top], '--c', linewidth=eq_linewidth)
    if detrended:
        ax3 = axarr[2].twinx()
        ax3.plot(detrended.dtarray, detrended.dU, marker='D', markersize=1.0, color='red', linestyle='')
        ax3.set_ylabel(detrended_label+' (mm)', fontsize=label_fontsize - 2, color='red')
        ax3.tick_params(labelcolor='red', labelsize=label_fontsize, axis='both')
    axarr[2].set_xlim([min(ts_obj.dtarray)-dt.timedelta(days=buffer_days),
                       max(ts_obj.dtarray) + dt.timedelta(days=buffer_days)])
    axarr[2].tick_params(labelsize=label_fontsize, rotation=label_rotation)

    axarr[0].set_title(title, fontsize=label_fontsize + 2)
    plt.savefig(savename, dpi=dpival, bbox_inches='tight')
    print("Saving figure as %s " % savename)
    return


def get_figure_name(plot_params, db_params):
    """
    Building filename and title strings: Metadata like station, offsets/outliers/seasonals, Datasource and Refframe.
    Title and Filename are programmatically linked.
    """
    if plot_params is None and db_params is None:
        raise ValueError("Error! Some kind of title and filename is required, "
                         "it cannot be automatically generated from nothing.")

    savename = db_params["station"]
    title = db_params["station"]

    title = title + ', ' + db_params["datasource"] + ' ' + db_params["refframe"]
    if plot_params["earthquakes_remove"] == 0 and plot_params["offsets_remove"] == 0 and \
            plot_params["seasonals_remove"] == 0:
        title = title + ', unaltered'
    if plot_params["earthquakes_remove"]:
        savename = savename + "_noeq"
        title = title + ', no earthquakes'
    if plot_params["seasonals_remove"]:  # If we are removing seasonals:
        savename = savename + "_noseasons"
        title = title + ', no seasonals'
        if plot_params["seasonals_type"] == "lssq":
            savename = savename + "_lssq"
            title = title + ' by least squares'
        elif plot_params["seasonals_type"] == "notch":
            savename = savename + "_notch"
            title = title + ' by notch filter'
        elif plot_params["seasonals_type"] == "grace":
            savename = savename + "_grace"
            title = title + ' by GRACE model'
        elif plot_params["seasonals_type"] == "stl":
            savename = savename + "_stl"
            title = title + ' by STL'
        elif plot_params["seasonals_type"] == "nldas":
            savename = savename + "_nldas"
            title = title + ' by NLDAS'
        elif plot_params["seasonals_type"] == "gldas":
            savename = savename + "_gldas"
            title = title + ' by GLDAS'
        elif plot_params["seasonals_type"] == "lsdm":
            savename = savename + "_lsdm"
            title = title + ' by LSDM'
        elif plot_params["seasonals_type"] == "shasta":
            savename = savename + "_shasta"
            title = title + ' by Shasta'
        elif plot_params["seasonals_type"] == "oroville":
            savename = savename + "_oroville"
            title = title + ' by Oroville'
        else:
            print("Error! Type of seasonal removal not recognized.")
    savename = savename + '_' + db_params["datasource"] + '_' + db_params["refframe"]
    savename = savename + "_ts.png"
    return title, savename
