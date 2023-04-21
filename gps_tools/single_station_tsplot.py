"""
Make a basic python plot of single-station position time series with corrections optional.
"""

import matplotlib.pyplot as plt
import subprocess
from . import gps_seasonal_removals, offsets, load_gnss


def view_single_station(station_name, data_config_file, offsets_remove=1, earthquakes_remove=0, outliers_remove=0,
                        seasonals_remove=0, starttime=None, endtime=None, outliers_def=15,
                        seasonals_type='lssq', datasource='cwu', refframe='NA', outdir=""):
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
    :param refframe: string, default NA
    :param outdir: string
    """
    db_params = config_db(station_name, datasource, refframe);
    plot_params = config_view(offsets_remove, earthquakes_remove, outliers_remove, outliers_def, seasonals_remove,
                              seasonals_type, starttime, endtime);
    [myData, offset_obj, eq_obj] = input_data(station_name, data_config_file, db_params);
    [updatedData, detrended] = compute(data_config_file, myData, offset_obj, eq_obj, plot_params);
    single_ts_plot(updatedData, detrended, plot_params, db_params, outdir);


# -------------- CONFIGURE ------------ #
def config_db(station, datasource, refframe):
    database_params = {"station": station, 'datasource': datasource, 'refframe': refframe};
    print("\n---------- VIEWING GNSS TIMESERIES ------------");
    print("Database Parameters: station = %s, dattasource = %s, refframe = %s " % (station, datasource, refframe));
    return database_params;


def config_view(offsets_remove, earthquakes_remove, outliers_remove, outliers_def, seasonals_remove, seasonals_type,
                starttime, endtime):
    """
    outliers_def : mm away from median filter.
    offsets_remove, earthquakes_remove, outliers_remove, seasonals_remove : booleans
    seasonals type : lssq, nldas, gldas, grace, lsdm, shasta, stl, notch
    """
    view_params = locals();  # package all provided arguments into a dictionary
    print("Viewing Parameters: earthquakes_remove = %d, outliers_remove = %d, seasonals_remove = %d"
          % (earthquakes_remove, outliers_remove, seasonals_remove));
    print("----------------------------------");
    return view_params;


# ----------- INPUTS ---------------- # 
def input_data(st_name, data_config_file, db_params):
    database = load_gnss.create_station_repo(data_config_file, db_params['refframe'], db_params['datasource']);
    [myData, offset_obj, eq_obj] = database.load_station(st_name);
    myData = myData.embed_tsobject_with_eqdates(eq_obj);  # embed data with eq object metadata
    return [myData, offset_obj, eq_obj];


# -------------- COMPUTE ------------ # 
def compute(data_config_file, myData, offset_obj, eq_obj, plot_params):
    if plot_params["starttime"] is None:
        plot_params["starttime"] = myData.dtarray[0];
    if plot_params["endtime"] is None:
        plot_params["endtime"] = myData.dtarray[-1];
    starttime = plot_params["starttime"]
    endtime = plot_params["endtime"]
    newData = myData.impose_time_limits(starttime, endtime);
    if plot_params["offsets_remove"]:  # Remove offsets and antenna changes
        newData = offsets.remove_offsets(newData, offset_obj);
    if plot_params["outliers_remove"]:  # Remove outliers
        outliers_def = plot_params["outliers_def"]
        newData = newData.remove_outliers(outliers_def);
    if plot_params["earthquakes_remove"]:  # Remove earthquakes
        newData = offsets.remove_offsets(newData, eq_obj);
    trend_out = gps_seasonal_removals.make_detrended_ts(newData, plot_params["seasonals_remove"],
                                                        plot_params["seasonals_type"], data_config_file);
    return [newData, trend_out];


# -------------- OUTPUTS ------------ # 
def single_ts_plot(ts_obj, detrended, plot_params, db_params, outdir):
    title, savename = get_figure_name(plot_params, db_params, outdir);
    label_fontsize = 18;

    # The major figure
    dpival = 500;
    # noinspection PyTypeChecker
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

    axarr[0].set_title(title, fontsize=label_fontsize + 2);
    plt.savefig(savename, dpi=dpival);
    print("Saving figure as %s " % savename)
    return;


def get_figure_name(plot_params, db_params, outdir):
    """
    Things that might go into the filename and title: Station, Offsets/outliers/seasonals, Datasource and Refframe
    """

    savename = outdir + db_params["station"];
    title = db_params["station"];
    if outdir != '':
        subprocess.call(['mkdir', '-p', outdir], shell=False);

    title = title + ', ' + db_params["datasource"] + ' ' + db_params["refframe"];
    if plot_params["earthquakes_remove"] == 0 and plot_params["offsets_remove"] == 0 and \
            plot_params["seasonals_remove"] == 0:
        title = title + ', unaltered';
    if plot_params["earthquakes_remove"]:
        savename = savename + "_noeq";
        title = title + ', no earthquakes'
    if plot_params["seasonals_remove"]:  # If we are removing seasonals:
        savename = savename + "_noseasons";
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
            savename = savename + "_nldas";
            title = title + ' by NLDAS'
        elif plot_params["seasonals_type"] == "gldas":
            savename = savename + "_gldas";
            title = title + ' by GLDAS'
        elif plot_params["seasonals_type"] == "lsdm":
            savename = savename + "_lsdm";
            title = title + ' by LSDM'
        elif plot_params["seasonals_type"] == "shasta":
            savename = savename + "_shasta";
            title = title + ' by Shasta';
        elif plot_params["seasonals_type"] == "oroville":
            savename = savename + "_oroville";
            title = title + ' by Oroville';
        else:
            print("Error! Type of seasonal removal not recognized.");
    savename = savename + '_' + db_params["datasource"] + '_' + db_params["refframe"];
    savename = savename + "_ts.png";
    return title, savename;
