"""
Plotting tools for GNSS stacks and movies
Plot stacks with trends, with steps, or anything else.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib, collections, scipy.ndimage
import matplotlib.cm as cm
import datetime as dt

StackParams = collections.namedtuple("StackParams", ['expname', 'bbox', 'eqtimes', 'starttime',
                                                     'endtime', 'labeltime', 'outdir', 'outname']);


def horizontal_full_ts(dataobj_list, distances, myparams, label="", removemean=1,
                       start_time_plot=dt.datetime.strptime("20050101", "%Y%m%d"),
                       end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d"),
                       label_date=dt.datetime.strptime("20200215", "%Y%m%d")):
    [_f, axarr] = plt.subplots(1, 2, sharex='all', sharey='all', figsize=(10, 8), dpi=300)

    spacing = 15;
    labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname, distances);

    color_boundary_object = matplotlib.colors.Normalize(vmin=closest_station, vmax=farthest_station, clip=True);
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
        _l1 = axarr[0].plot_date(dataobj_list[i].dtarray, edata, marker='+', markersize=1.5, color=line_color);
        axarr[0].text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    axarr[0].set_xlim(start_time_plot, end_time_plot);
    bottom, top = axarr[0].get_ylim();
    for i in range(len(myparams.eqtimes)):
        axarr[0].plot_date([myparams.eqtimes[i], myparams.eqtimes[i]], [bottom, top], '--k', linewidth=0.75);
    for i in range(len(labeltimes)):
        axarr[0].text(labeltimes[i], top - spacing / 2, labels[i], fontsize=10, color='blue', fontweight='bold');
    axarr[0].tick_params(axis='x', labelrotation=40)
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
        _l1 = axarr[1].plot_date(dataobj_list[i].dtarray, ndata, marker='+', markersize=1.5, color=line_color);
    axarr[1].set_xlim(start_time_plot, end_time_plot);
    bottom, top = axarr[1].get_ylim();
    for i in range(len(myparams.eqtimes)):
        axarr[1].plot_date([myparams.eqtimes[i], myparams.eqtimes[i]], [bottom, top], '--k', linewidth=0.75);
    axarr[1].set_ylabel("North (mm)");
    axarr[1].set_title("North GPS Time Series")
    axarr[1].grid(True)
    axarr[1].tick_params(axis='x', labelrotation=40)

    for i in range(len(labeltimes)):
        axarr[1].text(labeltimes[i], top - spacing / 2, labels[i], fontsize=10, color='blue', fontweight='bold');
    custom_cmap.set_array(range(int(closest_station), int(farthest_station)))
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + '.png')
    plt.close();
    print("Horizontal plots created.");

    return;


def vertical_full_ts(dataobj_list, distances, myparams, label="", removemean=1,
                     start_time_plot=dt.datetime.strptime("20050101", "%Y%m%d"),
                     end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d"),
                     label_date=dt.datetime.strptime("20200215", "%Y%m%d")):
    plt.figure(figsize=(6, 8), dpi=160);

    spacing = 40;
    labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname, distances);
    color_boundary_object = matplotlib.colors.Normalize(vmin=closest_station, vmax=farthest_station, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');

    # Vertical
    for i in range(len(dataobj_list)):
        offset = spacing * i;
        udata = dataobj_list[i].dU;
        umean = np.nanmean(dataobj_list[i].dU)
        if removemean == 0:
            umean = 0;  # leave the data un-meaned
        udata = [x + offset - umean for x in udata];
        line_color = custom_cmap.to_rgba(distances[i]);
        _l1 = plt.gca().plot_date(dataobj_list[i].dtarray, udata, marker='+', markersize=1.5, color=line_color);
        plt.gca().text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    plt.gca().set_xlim(start_time_plot, end_time_plot);
    plt.gca().tick_params(axis='x', labelrotation=40)
    bottom, top = plt.gca().get_ylim();
    for i in range(len(myparams.eqtimes)):
        plt.gca().plot_date([myparams.eqtimes[i], myparams.eqtimes[i]], [bottom, top], '--k', linewidth=0.75);
    plt.gca().set_ylabel("Vertical (mm)");
    plt.gca().set_title("Vertical GPS Time Series")
    plt.gca().grid(True)

    custom_cmap.set_array(range(int(closest_station), int(farthest_station)));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + '_vertical.png');
    plt.close();
    print("Vertical plot created.");
    return;


def horizontal_filtered_plots(dataobj_list, distances, myparams, label="",
                              start_time_plot=dt.datetime.strptime("20050101", "%Y%m%d"),
                              end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d"),
                              label_date=dt.datetime.strptime("20200215", "%Y%m%d")):
    [_f, axarr] = plt.subplots(2, 1, sharex='all', figsize=(15, 8), dpi=160);
    offset = -15;
    spacing = 2;

    labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname, distances);
    color_boundary_object = matplotlib.colors.Normalize(vmin=closest_station, vmax=farthest_station, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet_r');

    # East and North
    for i in range(len(dataobj_list)):
        offset = offset + spacing;
        umean = dataobj_list[i].dE[0];  # start at the beginning
        line_color = custom_cmap.to_rgba(distances[i]);
        # l1 = axarr[0].plot(dataobj_list[i].dtarray, dataobj_list[i].dU-umean, linestyle='solid', linewidth=0,
        #                    marker='.', color='red'); # debugging
        udata = scipy.ndimage.median_filter(dataobj_list[i].dE - umean, size=365);
        axarr[0].plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
        axarr[0].text(label_date, offset, dataobj_list[i].name, fontsize=9, color=line_color);
    axarr[0].set_xlim(start_time_plot, end_time_plot);
    bottom, top = axarr[0].get_ylim();
    for i in range(len(myparams.eqtimes)):
        axarr[0].plot_date([myparams.eqtimes[i], myparams.eqtimes[i]], [bottom, top], '--k');
    axarr[0].set_ylabel("Filtered East (mm)");
    axarr[0].set_title("Filtered GPS Time Series")
    axarr[0].grid(True)

    for i in range(len(dataobj_list)):
        umean = dataobj_list[i].dN[0];  # start at the beginning
        line_color = custom_cmap.to_rgba(distances[i]);
        # l1 = axarr[1].plot(dataobj_list[i].dtarray, dataobj_list[i].dU-umean, linestyle='solid', linewidth=0,
        #                    marker='.', color='red'); # debugging
        udata = scipy.ndimage.median_filter(dataobj_list[i].dN - umean, size=365);
        _l1 = axarr[1].plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
    bottom, top = axarr[1].get_ylim();
    for i in range(len(myparams.eqtimes)):
        axarr[1].plot_date([myparams.eqtimes[i], myparams.eqtimes[i]], [bottom, top], '--k');
    axarr[1].set_ylabel("Filtered North (mm)");
    axarr[1].grid(True)

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + '_horiz_filt.png');
    plt.close();
    print("Horiz plot created.");
    return;


def vertical_filtered_plots(dataobj_list, distances, myparams, label="",
                            start_time_plot=dt.datetime.strptime("20050101", "%Y%m%d"),
                            end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d"),
                            label_date=dt.datetime.strptime("20200215", "%Y%m%d")):
    plt.figure(figsize=(15, 8), dpi=160);

    labeltimes, labels, closest_station, farthest_station = configure_beautiful_plots(myparams.expname, distances);
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
        _l1 = plt.gca().plot(dataobj_list[i].dtarray, udata, linestyle='solid', linewidth=1, color=line_color);
        plt.gca().text(label_date, label_mm, dataobj_list[i].name, fontsize=9, color=line_color);
        offset = offset + 1.6;
    plt.gca().set_xlim(start_time_plot, end_time_plot);
    bottom, top = plt.gca().get_ylim();
    for i in range(len(myparams.eqtimes)):
        plt.gca().plot_date([myparams.eqtimes[i], myparams.eqtimes[i]], [bottom, top], '--k');
    plt.gca().set_ylabel("Filtered Vertical (mm)");
    plt.gca().set_title("Filtered Vertical GPS Time Series")
    plt.gca().grid(True)

    custom_cmap.set_array(range(int(closest_station), int(farthest_station)));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilometers from center');

    plt.savefig(myparams.outdir + "/" + myparams.outname + '_TS_' + label + 'vert_filt.png');
    plt.close();
    print("Vertical plot created.");
    return;


def configure_beautiful_plots(expname, distances):
    """
    Some hard-coded information to make beautiful plots for specific experiments.  Reminder for old code:

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
    """
    labeltimes, labels = [], [];

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

    return labeltimes, labels, closest_station, farthest_station;


def pygmt_map(dataobj_list, center, outdir, outname):
    """
    Optionally, use pygmt library to make a plot of stations locations
    """
    import pygmt
    offset = 0.2;

    lons = [x.coords[0] for x in dataobj_list];
    lats = [x.coords[1] for x in dataobj_list];
    names = [x.name for x in dataobj_list];
    region = [min(lons) - offset, max(lons) + offset, min(lats) - offset, max(lats) + offset];

    fig = pygmt.Figure()
    fig.basemap(region=region, projection="M8i", frame="0.25");
    # fig.grdimage("@earth_relief_30s",region=region,I="+d");  # takes a while the first time, but faster afterwards
    fig.coast(shorelines="0.5p,black", land='peachpuff2', water='skyblue', resolution="h");
    fig.coast(borders='1', shorelines='1.0p,black');
    fig.coast(borders='2', shorelines='0.5p,black');
    fig.text(x=lons, y=lats, text=names, font='15p,Helvetica-Bold,black', offset="0.24i/0.11i");
    fig.plot(x=lons, y=lats, style='c0.1i', color='black', pen='0.5p,black')
    fig.plot(x=center[0], y=center[1], style='a0.1i', color='red', pen='0.5p,red')
    fig.savefig(outdir + "/" + outname + '_map.png');
    print("Saving map %s" % (outdir + "/" + outname + '_map.png') );
    return;


def write_stack_params(myparams):
    outfilename = myparams.outdir+"/"+myparams.outname+"_params.txt";
    print('Writing param file %s ' % outfilename);
    ofile = open(outfilename, 'w');
    for name, value in zip(myparams._fields, myparams):
        ofile.write(name);
        ofile.write(": ");
        ofile.write(str(value));
        ofile.write("\n");
    ofile.close();
    return;
