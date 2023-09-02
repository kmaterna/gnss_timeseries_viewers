"""
Plotting tools for GNSS stacks and movies
Plot stacks with trends, with steps, or anything else.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib, scipy
import matplotlib.cm as cm
import datetime as dt
import json
from .gps_ts_functions import get_starttime, get_endtime


def write_params(outfile, param_dict):
    print("Writing stack parameters in %s " % outfile);
    with open(outfile, 'w') as fp:
        json.dump(param_dict, fp, indent=4, default=str)
    return;


def horizontal_full_ts(dataobj_list, distances, outname, removemean=1, start_time_plot=None,
                       end_time_plot=None, label_date=None, vmin=None, vmax=None, spacing=15, eqtimes=()):

    [_f, axarr] = plt.subplots(1, 2, sharex='all', sharey='all', figsize=(10, 8), dpi=300)
    start_plot, end_plot, label_date, vmin, vmax = get_plot_params(dataobj_list, start_time_plot, end_time_plot,
                                                                   label_date, vmin, vmax, distances);

    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');

    # East
    for i in range(len(dataobj_list)):
        offset = spacing * i;
        edata = dataobj_list[i].dE;
        emean = np.nanmean(dataobj_list[i].dE);
        if removemean == 0:
            emean = 0;  # leave the data un-meaned
        edata = [x + offset - emean for x in edata];
        line_color = custom_cmap.to_rgba(distances[i]);
        _l1 = axarr[0].plot_date(dataobj_list[i].dtarray, edata, marker='+', markersize=4.5, color=line_color);
        axarr[0].text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    axarr[0].set_xlim(start_plot, end_plot);
    bottom, top = axarr[0].get_ylim();
    for i in range(len(eqtimes)):
        axarr[0].plot_date([eqtimes[i], eqtimes[i]], [bottom, top], '--k', linewidth=0.75);
    axarr[0].tick_params(axis='x', labelrotation=40)
    axarr[0].set_ylabel("East (mm)");
    axarr[0].set_title("East GNSS Time Series")
    axarr[0].grid(True)

    # North
    for i in range(len(dataobj_list)):
        offset = spacing * i;
        ndata = dataobj_list[i].dN;
        nmean = np.nanmean(dataobj_list[i].dN);
        if removemean == 0:
            nmean = 0;  # leave the data un-meaned
        ndata = [x + offset - nmean for x in ndata];
        line_color = custom_cmap.to_rgba(distances[i]);
        _l1 = axarr[1].plot_date(dataobj_list[i].dtarray, ndata, marker='+', markersize=4.5, color=line_color);
    axarr[1].set_xlim(start_plot, end_plot);
    bottom, top = axarr[1].get_ylim();
    for i in range(len(eqtimes)):
        axarr[1].plot_date([eqtimes[i], eqtimes[i]], [bottom, top], '--k', linewidth=0.75);
    axarr[1].set_ylabel("North (mm)");
    axarr[1].set_title("North GNSS Time Series")
    axarr[1].grid(True)
    axarr[1].tick_params(axis='x', labelrotation=40)

    custom_cmap.set_array(range(int(vmin), int(vmax)));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    plt.savefig(outname);
    plt.close();
    print("Horizontal plots created in file %s" % outname);
    return;


def vertical_full_ts(dataobj_list, distances, outname, removemean=1, start_time_plot=None,
                     end_time_plot=None, label_date=None, vmin=None, vmax=None, spacing=40, eqtimes=()):
    plt.figure(figsize=(6, 8), dpi=160);
    start_plot, end_plot, label_date, vmin, vmax = get_plot_params(dataobj_list, start_time_plot, end_time_plot,
                                                                   label_date, vmin, vmax, distances);

    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');

    for i in range(len(dataobj_list)):  # plotting vertical
        offset = spacing * i;
        udata = dataobj_list[i].dU;
        umean = np.nanmean(dataobj_list[i].dU)
        if removemean == 0:
            umean = 0;  # leave the data un-meaned
        udata = [x + offset - umean for x in udata];
        line_color = custom_cmap.to_rgba(distances[i]);
        _l1 = plt.gca().plot_date(dataobj_list[i].dtarray, udata, marker='+', markersize=1.5, color=line_color);
        plt.gca().text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    plt.gca().set_xlim(start_plot, end_plot);
    plt.gca().tick_params(axis='x', labelrotation=40)
    bottom, top = plt.gca().get_ylim();
    for i in range(len(eqtimes)):
        plt.gca().plot_date([eqtimes[i], eqtimes[i]], [bottom, top], '--k', linewidth=0.75);
    plt.gca().set_ylabel("Vertical (mm)");
    plt.gca().set_title("Vertical GNSS Time Series")
    plt.gca().grid(True)

    custom_cmap.set_array(range(int(vmin), int(vmax)));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    plt.savefig(outname);
    plt.close();
    print("Vertical plot created in %s." % outname);
    return;


def horizontal_filtered_plots(dataobj_list, distances, outname, start_time_plot=None, end_time_plot=None,
                              label_date=None, vmin=None, vmax=None, spacing=2, offset=-15, eqtimes=()):
    [_f, axarr] = plt.subplots(2, 1, sharex='all', figsize=(15, 8), dpi=160);
    start_plot, end_plot, label_date, vmin, vmax = get_plot_params(dataobj_list, start_time_plot, end_time_plot,
                                                                   label_date, vmin, vmax, distances);

    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet_r');

    # East and North
    for i in range(len(dataobj_list)):
        offset = offset + spacing;
        umean = dataobj_list[i].dE[0];  # start at the beginning
        line_color = custom_cmap.to_rgba(distances[i]);
        udata = scipy.ndimage.median_filter(dataobj_list[i].dE - umean, size=365);
        axarr[0].plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
        axarr[0].text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    axarr[0].set_xlim(start_plot, end_plot);
    bottom, top = axarr[0].get_ylim();
    for i in range(len(eqtimes)):
        axarr[0].plot_date([eqtimes[i], eqtimes[i]], [bottom, top], '--k');
    axarr[0].set_ylabel("Filtered East (mm)");
    axarr[0].set_title("Filtered GNSS Time Series")
    axarr[0].grid(True)

    for i in range(len(dataobj_list)):
        umean = dataobj_list[i].dN[0];  # start at the beginning
        line_color = custom_cmap.to_rgba(distances[i]);
        # l1 = axarr[1].plot(dataobj_list[i].dtarray, dataobj_list[i].dU-umean, linestyle='solid', linewidth=0,
        #                    marker='.', color='red'); # debugging
        udata = scipy.ndimage.median_filter(dataobj_list[i].dN - umean, size=365);
        _l1 = axarr[1].plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
    axarr[1].set_xlim(start_plot, end_plot);
    bottom, top = axarr[1].get_ylim();
    for i in range(len(eqtimes)):
        axarr[1].plot_date([eqtimes[i], eqtimes[i]], [bottom, top], '--k');
    axarr[1].set_ylabel("Filtered North (mm)");
    axarr[1].grid(True)

    plt.savefig(outname);
    plt.close();
    print("Horiz plot created in %s." % outname);
    return;


def vertical_filtered_plots(dataobj_list, distances, outname, start_time_plot=None, end_time_plot=None,
                            label_date=None, vmin=None, vmax=None, spacing=1.6, offset=-10, eqtimes=()):
    plt.figure(figsize=(15, 8), dpi=160);
    start_plot, end_plot, label_date, vmin, vmax = get_plot_params(dataobj_list, start_time_plot, end_time_plot,
                                                                   label_date, vmin, vmax, distances);

    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet_r');

    # Vertical
    for i in range(len(dataobj_list)):
        # umean = dataobj_list[i].dU[0];  # start at the beginning for the trends-in vertical plot  # if trends-in
        # label_mm = dataobj_list[i].dU[-1] - dataobj_list[i].dU[0];  # if trends-in
        # udata = scipy.ndimage.median_filter(dataobj_list[i].dU - umean, size=60);  # if trends-in
        umean = np.mean(dataobj_list[i].dU);  # start at the mean for the detrended vertical plot
        label_mm = offset;
        udata = scipy.ndimage.median_filter(dataobj_list[i].dU - umean, size=365);

        line_color = custom_cmap.to_rgba(distances[i]);
        _l1 = plt.gca().plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
        plt.gca().text(label_date, label_mm, dataobj_list[i].name, fontsize=9, color=line_color);
        offset = offset + spacing;
    plt.gca().set_xlim(start_plot, end_plot);
    bottom, top = plt.gca().get_ylim();
    for i in range(len(eqtimes)):
        plt.gca().plot_date([eqtimes[i], eqtimes[i]], [bottom, top], '--k');
    plt.gca().set_ylabel("Filtered Vertical (mm)");
    plt.gca().set_title("Filtered Vertical GNSS Time Series")
    plt.gca().grid(True)

    custom_cmap.set_array(range(int(vmin), int(vmax)));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    plt.savefig(outname);
    plt.close();
    print("Vertical plot created in %s ." % outname);
    return;


def get_plot_params(ts_objects_list, start_time_plot, end_time_plot, label_date, vmin, vmax, color_array):
    """Override default parameters for plot start, plot end, and plot label, and plot colors"""
    if len(ts_objects_list) == 0:
        raise ValueError("Error! Length of timeseries objects provided to stacking plot is zero!");
    if start_time_plot is None:
        start_time_plot = get_starttime(ts_objects_list);
    if end_time_plot is None:
        end_time_plot = get_endtime(ts_objects_list);
    if label_date is None:
        label_date = get_labeltime(start_time_plot, end_time_plot, fraction=0.025);
    if vmin is None:
        vmin = np.min(color_array);
    if vmax is None:
        vmax = np.max(color_array);
    if start_time_plot > end_time_plot:
        raise ValueError("Error! Start time of plot is after end time. You likely have an error.");
    return start_time_plot, end_time_plot, label_date, vmin, vmax;


def get_labeltime(starttime, endtime, fraction=1/10):
    """
    :param starttime: dt.datetime object
    :param endtime: dt.datetime object
    :param fraction: default 1/10
    :return: interval + 1/10 of the interval
    """
    interval = endtime - starttime
    labeltime = endtime + dt.timedelta(days=int(interval.days*fraction));
    return labeltime;
