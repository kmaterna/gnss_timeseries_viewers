# Plotting tools for stacks and movies
# Whether we are plotting stacks with trends, with steps, or anything else, 
# This is the common set of plotting tools we want to use. 
# ------------------------------- # 
# Contains:
# horizontal_full_ts(dataobj_list, distances, myparams, label)
# vertical_full_ts(dataobj_list, distances, myparams, label)
# horizontal_filtered_plots(dataobj_list, distances, myparams, label)
# vertical_filtered_plots(dataobj_list, distances, myparams, label)
# pygmt_map(dataobj_list, myparams)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import datetime as dt
import scipy.ndimage


def horizontal_full_ts(dataobj_list, distances, myparams, label="", removemean=1):
    fig = plt.figure(figsize=(20, 15), dpi=160);
    [f, axarr] = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 8))
    label_date = dt.datetime.strptime("20200215", "%Y%m%d");
    start_time_plot = dt.datetime.strptime("20050101", "%Y%m%d");
    end_time_plot = dt.datetime.strptime("20200116", "%Y%m%d");

    spacing = 15;
    EQtimes, labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname,
                                                                                               distances);

    color_boundary_object = matplotlib.colors.Normalize(vmin=closest_station, vmax=farthest_station, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet_r');

    # East
    for i in range(len(dataobj_list)):
        offset = spacing * i;
        edata = dataobj_list[i].dE;
        emean = np.nanmean(dataobj_list[i].dE);
        if removemean == 0:
            emean = 0;  # leave the data un-meaned
        edata = [x + offset - emean for x in edata];
        line_color = custom_cmap.to_rgba(distances[i]);
        l1 = axarr[0].plot_date(dataobj_list[i].dtarray, edata, marker='+', markersize=1.5, color=line_color);
        axarr[0].text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    axarr[0].set_xlim(start_time_plot, end_time_plot);
    # axarr[0].set_ylim([-20,offset+15])
    bottom, top = axarr[0].get_ylim();
    for i in range(len(EQtimes)):
        axarr[0].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k');
    for i in range(len(labeltimes)):
        axarr[0].text(labeltimes[i], top - spacing / 2, labels[i], fontsize=10, color='blue', fontweight='bold');
    axarr[0].set_ylabel("East (mm)");
    axarr[0].set_title("East GPS Time Series")
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
        l1 = axarr[1].plot_date(dataobj_list[i].dtarray, ndata, marker='+', markersize=1.5, color=line_color);
    axarr[1].set_xlim(start_time_plot, end_time_plot);
    # axarr[1].set_ylim([-20,offset+15])
    bottom, top = axarr[1].get_ylim();
    for i in range(len(EQtimes)):
        axarr[1].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k');
    axarr[1].set_ylabel("North (mm)");
    axarr[1].set_title("North GPS Time Series")
    axarr[1].grid(True)
    for i in range(len(labeltimes)):
        axarr[1].text(labeltimes[i], top - spacing / 2, labels[i], fontsize=10, color='blue', fontweight='bold');
    custom_cmap.set_array(range(int(closest_station), int(farthest_station)))
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + '.png')
    plt.close();
    print("Horizontal plots created.");

    return;


def vertical_full_ts(dataobj_list, distances, myparams, label="", removemean=1):
    plt.figure(figsize=(6, 8), dpi=160);
    label_date = dt.datetime.strptime("20200215", "%Y%m%d");
    start_time_plot = dt.datetime.strptime("20050101", "%Y%m%d");
    end_time_plot = dt.datetime.strptime("20200116", "%Y%m%d");

    spacing = 40;
    EQtimes, labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname,
                                                                                               distances);
    color_boundary_object = matplotlib.colors.Normalize(vmin=closest_station, vmax=farthest_station, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet_r');

    # Vertical
    for i in range(len(dataobj_list)):
        offset = spacing * i;
        udata = dataobj_list[i].dU;
        umean = np.nanmean(dataobj_list[i].dU)
        if removemean == 0:
            umean = 0;  # leave the data un-meaned
        udata = [x + offset - umean for x in udata];
        line_color = custom_cmap.to_rgba(distances[i]);
        l1 = plt.gca().plot_date(dataobj_list[i].dtarray, udata, marker='+', markersize=1.5, color=line_color);
        plt.gca().text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    plt.gca().set_xlim(start_time_plot, end_time_plot);
    # plt.gca().set_ylim([-20,offset+20])
    bottom, top = plt.gca().get_ylim();
    for i in range(len(EQtimes)):
        plt.gca().plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k');
    plt.gca().set_ylabel("Vertical (mm)");
    plt.gca().set_title("Vertical GPS Time Series")
    plt.gca().grid(True)

    custom_cmap.set_array(range(int(closest_station), int(farthest_station)));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    # # new axis for plotting the map of california
    # ax = plt.axes(position=[0.8, 0.1, 0.2, 0.2], xticklabels=[], yticklabels=[]);
    # [ca_lons, ca_lats] = np.loadtxt('../california_bdr', unpack=True);
    # ax.plot(ca_lons, ca_lats, 'k');
    # for i in range(len(dataobj_list)):
    #     ax.plot(dataobj_list[i].coords[0], dataobj_list[i].coords[1], '.g', markersize=0.6);

    # new axis for extra labels
    ax = plt.axes(position=[0.8, 0.4, 0.2, 0.2], xticklabels=[], yticklabels=[]);
    ax.set_ylim([-1, 1])
    ax.set_xlim([-0.1, 1])
    ax.text(0, 0.75, myparams.proc_center + " data");
    ax.text(0, 0.37, myparams.center);
    ax.text(0, 0, str(myparams.radius) + " km radius");

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + '_vertical.png');
    plt.close();
    print("Vertical plot created.");
    return;


def horizontal_filtered_plots(dataobj_list, distances, myparams, label=""):
    [f, axarr] = plt.subplots(2, 1, sharex=True, figsize=(15, 8), dpi=160);
    label_date = dt.datetime.strptime("20200215", "%Y%m%d");
    start_time_plot = dt.datetime.strptime("20050101", "%Y%m%d");
    end_time_plot = dt.datetime.strptime("20200116", "%Y%m%d");
    offset = -15;
    spacing = 2;

    EQtimes, labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname,
                                                                                               distances);
    color_boundary_object = matplotlib.colors.Normalize(vmin=closest_station, vmax=farthest_station, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet_r');

    # East and North
    for i in range(len(dataobj_list)):
        offset = offset + spacing;
        # umean=np.mean(dataobj_list[i].dE);  # start at the mean.
        umean = dataobj_list[i].dE[0];  # start at the beginning
        line_color = custom_cmap.to_rgba(distances[i]);
        # l1 = plt.gca().plot(dataobj_list[i].dtarray,dataobj_list[i].dU-umean,linestyle='solid',linewidth=0,marker='.',color='red' ); # for debugging the filter
        udata = scipy.ndimage.median_filter(dataobj_list[i].dE - umean, size=365);
        axarr[0].plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
        axarr[0].text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    axarr[0].set_xlim(start_time_plot, end_time_plot);
    bottom, top = axarr[0].get_ylim();
    for i in range(len(EQtimes)):
        axarr[0].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k');
    axarr[0].set_ylabel("Filtered East (mm)");
    axarr[0].set_title("Filtered GPS Time Series")
    axarr[0].grid(True)

    for i in range(len(dataobj_list)):
        # umean=np.mean(dataobj_list[i].dE);  # start at the mean.
        umean = dataobj_list[i].dN[0];  # start at the beginning
        line_color = custom_cmap.to_rgba(distances[i]);
        # l1 = plt.gca().plot(dataobj_list[i].dtarray,dataobj_list[i].dU-umean,linestyle='solid',linewidth=0,marker='.',color='red' ); # for debugging the filter
        udata = scipy.ndimage.median_filter(dataobj_list[i].dN - umean, size=365);
        l1 = axarr[1].plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
    bottom, top = axarr[1].get_ylim();
    for i in range(len(EQtimes)):
        axarr[1].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k');
    axarr[1].set_ylabel("Filtered North (mm)");
    axarr[1].grid(True)

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + '_horiz_filt.png');
    plt.close();
    print("Horiz plot created.");
    return;


def vertical_filtered_plots(dataobj_list, distances, myparams, label=""):
    plt.figure(figsize=(15, 8), dpi=160);
    label_date = dt.datetime.strptime("20200215", "%Y%m%d");
    start_time_plot = dt.datetime.strptime("20050101", "%Y%m%d");
    end_time_plot = dt.datetime.strptime("20200116", "%Y%m%d");

    EQtimes, labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname,
                                                                                               distances);
    color_boundary_object = matplotlib.colors.Normalize(vmin=closest_station, vmax=farthest_station, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet_r');
    offset = -10;

    # Vertical
    for i in range(len(dataobj_list)):
        if label == "trendsin_":
            umean = dataobj_list[i].dU[0];  # start at the beginning for the trends-in vertical plot
            label_mm = dataobj_list[i].dU[-1] - dataobj_list[i].dU[0];
            udata = scipy.ndimage.median_filter(dataobj_list[i].dU - umean, size=60);
        else:
            umean = np.mean(dataobj_list[i].dU);  # start at the mean for the detrended vertical plot
            label_mm = offset;
            udata = scipy.ndimage.median_filter(dataobj_list[i].dU - umean, size=365);

        line_color = custom_cmap.to_rgba(distances[i]);
        l1 = plt.gca().plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
        plt.gca().text(label_date, label_mm, dataobj_list[i].name, fontsize=9, color=line_color);
        offset = offset + 1.6;
    plt.gca().set_xlim(start_time_plot, end_time_plot);
    bottom, top = plt.gca().get_ylim();
    for i in range(len(EQtimes)):
        plt.gca().plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k');
    plt.gca().set_ylabel("Filtered Vertical (mm)");
    plt.gca().set_title("Filtered Vertical GPS Time Series")
    plt.gca().grid(True)

    custom_cmap.set_array(range(int(closest_station), int(farthest_station)));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    # # new axis for plotting the map of california
    # ax = plt.axes(position=[0.8, 0.1, 0.1, 0.2], xticklabels=[], yticklabels=[]);
    # [ca_lons, ca_lats] = np.loadtxt('../california_bdr', unpack=True);
    # ax.plot(ca_lons, ca_lats, 'k');
    # for i in range(len(dataobj_list)):
    #     ax.plot(dataobj_list[i].coords[0], dataobj_list[i].coords[1], '.g', markersize=0.6);

    # new axis for extra labels
    ax = plt.axes(position=[0.8, 0.4, 0.1, 0.2], xticklabels=[], yticklabels=[]);
    ax.set_ylim([-1, 1])
    ax.set_xlim([-0.1, 1])
    ax.text(0, 0.75, myparams.proc_center + " data");
    ax.text(0, 0.37, myparams.center);
    ax.text(0, 0, str(myparams.radius) + " km radius");

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + 'vert_filt.png');
    plt.close();
    print("Vertical plot created.");
    return;


def configure_beautiful_plots(expname, distances):
    # Some hard coded information to make beautiful plots for specific experiments
    EQtimes = [];
    labeltimes = [];
    labels = [];

    if expname == 'Mend' or expname == "Humboldt":
        # What black lines do you want added to the figure? This is good for Mendocino
        EQtimes.append(dt.datetime.strptime("20140310", "%Y%m%d"));  # starts with the most important one
        EQtimes.append(dt.datetime.strptime("20050615", "%Y%m%d"));  # other earthquakes added to the figure
        EQtimes.append(dt.datetime.strptime("20100110", "%Y%m%d"));
        EQtimes.append(dt.datetime.strptime("20161208", "%Y%m%d"));

    elif expname == "SSGF":
        # # This is good for SoCal
        EQtimes.append(dt.datetime.strptime("20100403", "%Y%m%d"));  # starts with the most important one
        EQtimes.append(dt.datetime.strptime("20050615", "%Y%m%d"));  # other earthquakes added to the figure
        EQtimes.append(dt.datetime.strptime("20120808", "%Y%m%d"));

    if expname == "Mend":
        closest_station = 70;
        farthest_station = 120;
    else:
        closest_station = min(distances);  # km from event
        farthest_station = max(distances);  # km from event

    # Labels
    if expname == "Mend" or expname == "Humboldt":
        labeltimes.append(dt.datetime.strptime('20070701', "%Y%m%d"));
        labeltimes.append(dt.datetime.strptime('20120101', "%Y%m%d"));
        labeltimes.append(dt.datetime.strptime('20150501', "%Y%m%d"));
        labeltimes.append(dt.datetime.strptime('20171201', "%Y%m%d"));
        labels.append("T1");
        labels.append("T2");
        labels.append("T3");
        labels.append("T4");

    return EQtimes, labeltimes, labels, closest_station, farthest_station;


def pygmt_map(dataobj_list, myparams):
    # Optionally, use the pygmt library to make a plot of the stations themselves.
    import pygmt
    offset = 0.2;

    lons = [];
    lats = [];
    names = [];
    for i in range(len(dataobj_list)):
        lons.append(dataobj_list[i].coords[0]);
        lats.append(dataobj_list[i].coords[1]);
        names.append(dataobj_list[i].name);
    region = [min(lons) - offset, max(lons) + offset, min(lats) - offset, max(lats) + offset];

    fig = pygmt.Figure()
    fig.basemap(region=region, projection="M8i", B="0.25");
    # fig.grdimage("@earth_relief_30s",region=region,I="+d");  # takes a little while the first time, but faster each time afterwards
    fig.coast(shorelines="0.5p,black", G='peachpuff2', S='skyblue', D="h");
    fig.coast(N='1', W='1.0p,black');
    fig.coast(N='2', W='0.5p,black');
    fig.text(x=[i + 0.06 for i in lons], y=lats, text=names, font='15p,Helvetica-Bold,black');
    fig.plot(x=lons, y=lats, S='c0.1i', G='black', W='0.5p,black')
    fig.plot(x=myparams.center[0], y=myparams.center[1], S='a0.1i', G='red', W='0.5p,red')
    fig.savefig(myparams.outdir + "/" + myparams.outname + '_map.png');
    return;
