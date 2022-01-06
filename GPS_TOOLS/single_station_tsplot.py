"""
Make a basic python plot of single-station position time series with corrections optional.
Types of seasonal options:
   lssq: fits seasonals and linear trend by least squares inversion.
  notch: removes the 1-year and 6-month components by notch filter.
  grace: uses GRACE loading model interpolated between monthly points where available, and linear inversion where not.
    stl: uses a pre-computed look-up table for STL time series.
"""

import matplotlib.pyplot as plt
import collections
from . import gps_ts_functions, gps_seasonal_removals, offsets, gps_io_functions, gps_input_pipeline

# Parameters for controlling plotting
Parameters = collections.namedtuple("Parameters", ['station', 'outliers_remove', 'outliers_def',
                                                   'earthquakes_remove', 'offsets_remove', 'seasonals_remove',
                                                   'seasonals_type', 'datasource', 'refframe', 'data_config_file']);


def view_single_station(station_name, offsets_remove=1, earthquakes_remove=0, outliers_remove=0, seasonals_remove=0,
                        starttime=None, endtime=None, outliers_def=15,
                        seasonals_type='lssq', datasource='cwu', refframe='NA',
                        data_config_file="data_config.txt", outdir="single_plots/"):
    """
    :param station_name: string
    :param offsets_remove: bool, default True
    :param earthquakes_remove: bool, default False
    :param outliers_remove: bool, default False
    :param seasonals_remove: bool, default False
    :param starttime: dt.datetime object, default None
    :param endtime: dt.datetime object, default None
    :param outliers_def: mm, default 15
    :param seasonals_type: string, default 'lssq'
    :param datasource: string, default cwu
    :param refframe: string, default NA
    :param data_config_file: string
    :param outdir: string
    """
    MyParams = configure(station_name, offsets_remove, earthquakes_remove, outliers_remove, outliers_def,
                         seasonals_remove, seasonals_type, datasource, refframe, data_config_file);
    [myData, offset_obj, eq_obj] = input_data(MyParams.station, MyParams.datasource, MyParams.refframe,
                                              MyParams.data_config_file);
    [updatedData, detrended] = compute(myData, offset_obj, eq_obj, MyParams, starttime, endtime);
    single_ts_plot(updatedData, detrended, MyParams, outdir);


# -------------- CONFIGURE ------------ # 
def configure(station, offsets_remove, earthquakes_remove, outliers_remove, outliers_def,
              seasonals_remove, seasonals_type, datasource, refframe, data_config_file):
    """
    outliers_def : mm away from median filter.
    refframe : itrf or nam
    offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove : booleans
    seasonals type : lssq, nldas, gldas, grace, lsdm, shasta, stl, notch
    datasource : unr, cwu, pbo, nmt, usgs
    """
    MyParams = Parameters(station=station, outliers_remove=outliers_remove, outliers_def=outliers_def,
                          earthquakes_remove=earthquakes_remove,
                          offsets_remove=offsets_remove, seasonals_remove=seasonals_remove,
                          seasonals_type=seasonals_type,
                          datasource=datasource, refframe=refframe, data_config_file=data_config_file);
    print("------- %s --------" % station);
    print(
        "Viewing station %s, earthquakes_remove=%d, outliers_remove=%d, seasonals_remove=%d, datasource=%s, refframe=%s"
        % (station, earthquakes_remove, outliers_remove, seasonals_remove, datasource, refframe));
    return MyParams;


# ----------- INPUTS ---------------- # 
def input_data(st_name, datasource, refframe, data_config_file):
    [myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(st_name, datasource, data_config_file, refframe);
    eqdates = [x.evdts for x in eq_obj];
    # First, we embed the data with the eq object metadata (always useful)
    myData = gps_io_functions.Timeseries(name=myData.name, coords=myData.coords, dtarray=myData.dtarray, dN=myData.dN,
                                         dE=myData.dE, dU=myData.dU, Sn=myData.Sn, Se=myData.Se, Su=myData.Su,
                                         EQtimes=eqdates);
    return [myData, offset_obj, eq_obj];


# -------------- COMPUTE ------------ # 
def compute(myData, offset_obj, eq_obj, MyParams, starttime, endtime):
    if starttime is None:
        starttime = myData.dtarray[0];
    if endtime is None:
        endtime = myData.dtarray[-1];
    newData = gps_ts_functions.impose_time_limits(myData, starttime, endtime);
    if MyParams.offsets_remove == 1:  # Remove offsets and antenna changes
        newData = offsets.remove_offsets(newData, offset_obj);
    if MyParams.outliers_remove == 1:  # Remove outliers
        newData = gps_ts_functions.remove_outliers(newData, MyParams.outliers_def);
    if MyParams.earthquakes_remove == 1:  # Remove earthquakes
        newData = offsets.remove_offsets(newData, eq_obj);
    trend_out = gps_seasonal_removals.make_detrended_ts(newData, MyParams.seasonals_remove, MyParams.seasonals_type,
                                                        MyParams.data_config_file);
    return [newData, trend_out];


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, detrended, MyParams, outdir):
    label_fontsize = 18;

    # The major figure
    dpival = 500;
    [_, axarr] = plt.subplots(3, 1, sharex=True, figsize=(10, 7), dpi=dpival);
    axarr[0].plot_date(ts_obj.dtarray, ts_obj.dE, color='blue', markeredgecolor='black', markersize=1.5);
    axarr[0].grid(linestyle='--', linewidth=0.5);
    axarr[0].set_ylabel('east (mm)', fontsize=label_fontsize);
    bottom, top = axarr[0].get_ylim();
    for i in range(len(ts_obj.EQtimes)):
        axarr[0].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k', linewidth=0.5);
    ax1 = axarr[0].twinx();
    ax1.plot_date(detrended.dtarray, detrended.dE, marker='D', markersize=1.0, color='red');
    ax1.set_ylabel('detrended (mm)', fontsize=label_fontsize - 2, color='red');
    ax1.tick_params(labelcolor='red', labelsize=label_fontsize, axis='both')
    axarr[0].tick_params(labelsize=label_fontsize);

    axarr[1].plot_date(ts_obj.dtarray, ts_obj.dN, color='blue', markeredgecolor='black', markersize=1.5);
    axarr[1].grid(linestyle='--', linewidth=0.5);
    axarr[1].set_ylabel('north (mm)', fontsize=label_fontsize);
    bottom, top = axarr[1].get_ylim();
    for i in range(len(ts_obj.EQtimes)):
        axarr[1].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k', linewidth=0.5);
    ax2 = axarr[1].twinx();
    ax2.plot_date(detrended.dtarray, detrended.dN, marker='D', markersize=1.0, color='red');
    ax2.set_ylabel('detrended (mm)', fontsize=label_fontsize - 2, color='red');
    ax2.tick_params(labelcolor='red', labelsize=label_fontsize, axis='both')
    axarr[1].tick_params(labelsize=label_fontsize);

    axarr[2].plot_date(ts_obj.dtarray, ts_obj.dU, color='blue', markeredgecolor='black', markersize=1.5);
    axarr[2].grid(linestyle='--', linewidth=0.5);
    axarr[2].set_ylabel('vertical (mm)', fontsize=label_fontsize);
    bottom, top = axarr[2].get_ylim();
    for i in range(len(ts_obj.EQtimes)):
        axarr[2].plot_date([ts_obj.EQtimes[i], ts_obj.EQtimes[i]], [bottom, top], '--k', linewidth=0.5);
    ax3 = axarr[2].twinx();
    ax3.plot_date(detrended.dtarray, detrended.dU, marker='D', markersize=1.0, color='red');
    ax3.set_ylabel('detrended (mm)', fontsize=label_fontsize - 2, color='red');
    axarr[2].set_xlim([min(ts_obj.dtarray), max(ts_obj.dtarray)]);
    ax3.tick_params(labelcolor='red', labelsize=label_fontsize, axis='both')
    axarr[2].tick_params(labelsize=label_fontsize);

    title, savename = get_figure_name(MyParams, outdir);
    axarr[0].set_title(title, fontsize=label_fontsize + 2);
    plt.savefig(savename, dpi=dpival);
    print("Saving figure as %s " % savename)
    return;


def get_figure_name(MyParams, outdir):
    """
    Things that might go into the name:
    1. Station
    2. Earthquakes removed
    3. Outliers removed
    4. Seasonals removed
    5. Datasource
    6. Refframe
    """

    savename = outdir + MyParams.station;
    title = MyParams.station;

    title = title + ', ' + MyParams.datasource + ' ' + MyParams.refframe
    if MyParams.earthquakes_remove == 0 and MyParams.offsets_remove == 0 and MyParams.seasonals_remove == 0:
        title = title + ', unaltered';
    if MyParams.earthquakes_remove:
        savename = savename + "_noeq";
        title = title + ', no earthquakes'
    if MyParams.seasonals_remove:  # If we are removing seasonals:
        savename = savename + "_noseasons";
        title = title + ', no seasonals'
        if MyParams.seasonals_type == "lssq":
            savename = savename + "_lssq"
            title = title + ' by least squares'
        if MyParams.seasonals_type == "notch":
            savename = savename + "_notch"
            title = title + ' by notch filter'
        if MyParams.seasonals_type == "grace":
            savename = savename + "_grace"
            title = title + ' by GRACE model'
        if MyParams.seasonals_type == "stl":
            savename = savename + "_stl"
            title = title + ' by STL'
        if MyParams.seasonals_type == "nldas":
            savename = savename + "_nldas";
            title = title + ' by NLDAS'
        if MyParams.seasonals_type == "gldas":
            savename = savename + "_gldas";
            title = title + ' by GLDAS'
        if MyParams.seasonals_type == "lsdm":
            savename = savename + "_lsdm";
            title = title + ' by LSDM'
        if MyParams.seasonals_type == "shasta":
            savename = savename + "_shasta";
            title = title + ' by Shasta';
        if MyParams.seasonals_type == "oroville":
            savename = savename + "_oroville";
            title = title + ' by Oroville';
    savename = savename + '_' + MyParams.datasource + '_' + MyParams.refframe;
    savename = savename + "_ts.jpg";
    return title, savename;
