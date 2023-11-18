"""
Download GNSS time series, velocities, and offsets from the USGS (https://earthquake.usgs.gov/monitoring/gps/)
One sub-network at a time
Call this script from the main USGS_Data directory on your local machine. 
Script by K. Materna, 2020
"""


import requests
import numpy as np
import pandas as pd
import datetime as dt
import os
import subprocess

# GLOBAL VARIABLES
base_url = 'https://earthquake.usgs.gov/monitoring/gps/'
networks = ['BasinAndRange',
            'BasinAndRange_SGPS',
            'CentralCalifornia',
            'CentralCalifornia_SGPS',
            'CentralUS',
            'DiabloCanyon_SPGPS',
            'EasternORWA_SGPS',
            'ECSZ_SGPS',
            'ElCentro',
            'Helens',
            'LongValley',
            'Mammoth_SGPS',
            'MtBaker',
            'NCalifornia_SGPS',
            'NorthEasternCal_SGPS',
            'Northridge',
            'Pacific_Northwest',
            'SFBayArea',
            'SFBayArea_SGPS',
            'Sisters',
            'Southern_California',
            'WindKetchFlat_SGPS',
            'YellowstoneContin',
            'Yellowstone_SPGPS'];
Networks_ignored = ['Alaska',
                    'Alaska_SGPS',
                    'CentralCalifornia_ITRF2014',
                    'JuanDeFuca',
                    'LKY2',
                    'Pacific']
vel_base_directory = os.path.join("Velocities", "")   # where will velocity files live on the local machine?
ts_base_directory = os.path.join("Time_Series", "")   # where will time series files live on the local machine?
offsets_directory = os.path.join("Offsets", "")       # where will the offsets live on the local machine?
usgs_cache_file = os.path.join("Metadata", "usgs_station_cache.txt")  # cache necessary for coords/start/end times


def download_usgs_velocity_tables(network, base_directory, outfile):
    """
    Go online and grab USGS velocity tables, in NA, ITRF2008, and Filtered
    One table for each network, as organized in the USGS database.
    """
    url = base_url + network + '/velocities';
    df_list = pd.read_html(requests.get(url).content)
    print("Reading %d tables on website %s " % (len(df_list), url) );
    if len(df_list) == 3:
        os.makedirs(base_directory, exist_ok=True);
        print("Assuming they follow the pattern NAM, ITRF2008, Filtered");

        # NORTH AMERICA
        nam_outfile = base_directory + "NAM_" + outfile;
        one_outfile = open(nam_outfile, 'w');
        one_outfile.write("# USGS Velocity File in NAM\n")
        one_outfile.write("# Downloaded from %s on %s\n" % (url, dt.datetime.now()));
        one_outfile.write("# ");
        one_outfile.close();
        df_list[0].to_csv(nam_outfile, index=False, sep=' ', mode='a');
        print("Writing %s " % nam_outfile);

        # ITRF2008
        itrf_outfile = base_directory + "ITRF_" + outfile;
        one_outfile = open(itrf_outfile, 'w');
        one_outfile.write("# USGS Velocity File in ITRF2008\n")
        one_outfile.write("# Downloaded from %s on %s\n" % (url, dt.datetime.now()));
        one_outfile.write("# ");
        one_outfile.close();
        df_list[1].to_csv(itrf_outfile, index=False, sep=' ', mode='a');
        print("Writing %s " % itrf_outfile);

        # REGIONALLY FILTERED
        filt_outfile = base_directory + "FILT_" + outfile;
        one_outfile = open(filt_outfile, 'w');
        one_outfile.write("# USGS Velocity File Regionally Filtered\n")
        one_outfile.write("# Downloaded from %s on %s\n" % (url, dt.datetime.now()));
        one_outfile.write("# ");
        one_outfile.close();
        df_list[2].to_csv(filt_outfile, index=False, sep=' ', mode='a');
        print("Writing %s " % filt_outfile);

    else:
        print("Skipping network %s because of unexpected number of reference frames (%d)." % (network, len(df_list)) )
    return;


def download_usgs_offset_tables(network, base_directory, outfile):
    url = base_url + network + '/offsets';
    try:
        df_list = pd.read_html(requests.get(url).content);
    except ValueError:
        print("No tables found for %s " % network);
        return;
    print("Reading %d tables on website %s " % (len(df_list), url) );
    if len(df_list) == 3:
        os.makedirs(base_directory, exist_ok=True);
        print("Assuming they follow the pattern NAM, ITRF2008, Filtered");

        # NORTH AMERICA
        nam_outfile = base_directory + "NA_" + outfile;
        one_outfile = open(nam_outfile, 'w');
        one_outfile.write("# USGS Offset File in NAM\n")
        one_outfile.write("# Downloaded from %s on %s\n" % (url, dt.datetime.now()));
        one_outfile.write("# ");
        one_outfile.close();
        df_list[0].to_csv(nam_outfile, index=False, sep=' ', mode='a');
        print("Writing %s " % nam_outfile);

        # ITRF2008
        itrf_outfile = base_directory + "ITRF_" + outfile;
        one_outfile = open(itrf_outfile, 'w');
        one_outfile.write("# USGS Offset File in ITRF2008\n")
        one_outfile.write("# Downloaded from %s on %s\n" % (url, dt.datetime.now()));
        one_outfile.write("# ");
        one_outfile.close();
        df_list[1].to_csv(itrf_outfile, index=False, sep=' ', mode='a');
        print("Writing %s " % itrf_outfile);

        # REGIONALLY FILTERED
        filt_outfile = base_directory + "FILT_" + outfile;
        one_outfile = open(filt_outfile, 'w');
        one_outfile.write("# USGS Offset File Regionally Filtered\n")
        one_outfile.write("# Downloaded from %s on %s\n" % (url, dt.datetime.now()));
        one_outfile.write("# ");
        one_outfile.close();
        df_list[2].to_csv(filt_outfile, index=False, sep=' ', mode='a');
        print("Writing %s " % filt_outfile);

    else:
        print("Skipping network %s because of unexpected number of reference frames." % network)
    return;


def download_usgs_ts_file(station, network, ts_base_directory):
    """
    Download some time series files from USGS in NA and ITRF2008
    They will be organized by network, as they are in the USGS database
    """
    directory_name = os.path.join(ts_base_directory + network, "");
    os.makedirs(directory_name, exist_ok=True);
    url_name = base_url + 'data/networks/' + network + '/' + station.lower() + '/nafixed/' + station.lower() + '.rneu'
    subprocess.call(["wget", url_name, '-O', directory_name+station.lower()+'_NAfixed.rneu'], shell=False);
    url_name = base_url + 'data/networks/' + network + '/' + station.lower() + '/itrf2008/' + station.lower() + '.rneu'
    subprocess.call(["wget", url_name, '-O', directory_name+station.lower()+'_ITRF2008.rneu'], shell=False);
    return;


def get_usgs_network_directory_from_velfile(infile):
    """
    Network parsing and plumbing
    We assume a parallel directory above the Velocity directory with a bunch of TS files
    From which we can derive startdate and enddate
    """
    usgs_network = os.path.split(infile)[1][0:-9];
    if usgs_network[0:4] == 'NAM_':
        usgs_network = usgs_network[4:];
    if usgs_network[0:5] == 'ITRF_':
        usgs_network = usgs_network[5:];
    head, _ = os.path.split(infile);
    head, _ = os.path.split(head);
    head, _ = os.path.split(head);
    usgs_directory = os.path.join(head, "");
    network_ts_directory = os.path.join(usgs_directory, 'Time_Series', usgs_network, "");
    total_vel_directory = os.path.join(usgs_directory, 'Velocities', usgs_network, "");
    return network_ts_directory, total_vel_directory;


def cache_usgs_start_endtimes():
    print("Writing network cache %s " % usgs_cache_file );
    ofile = open(usgs_cache_file, 'w');
    for network in networks:
        velfile = os.path.join(vel_base_directory + network, 'NAM_' + network + '_vels.txt');
        network_ts_directory, _ = get_usgs_network_directory_from_velfile(velfile);

        [names, lon, lat] = np.loadtxt(velfile, skiprows=3, usecols=(0, 1, 2), unpack=True, dtype={'names': (
                'name', 'lon', 'lat'), 'formats': ('U4', float, float)});
        # Populating the first_epoch and last_epoch with information from the associated time series directory.
        for i in range(len(names)):
            ts_filename = os.path.join(network_ts_directory, names[i] + '_NAfixed.rneu');
            dates = np.loadtxt(ts_filename, unpack=True, usecols=(0,),
                               dtype={'names': ('dtstrs',), 'formats': ('U8',)}, ndmin=1);
            if len(dates) == 0:   # for the rare empty file
                continue;
            first_epoch = dates[0][0];
            last_epoch = dates[-1][0];
            ofile.write("%s %.5f %.5f %s %s %s\n" % (names[i], lon[i], lat[i], first_epoch, last_epoch, network) );

    ofile.close();
    return;


# THE SUB-DRIVERS
def download_usgs_offsets():
    for network in networks:
        directory_name = os.path.join(offsets_directory + network, '');
        download_usgs_offset_tables(network, directory_name, outfile=network+'_offsets.txt');
    return;


def download_usgs_velocities():
    for network in networks:
        directory_name = os.path.join(vel_base_directory + network, '');
        download_usgs_velocity_tables(network, directory_name, outfile=network+"_vels.txt");    
    return;


def download_usgs_time_series():
    for network in networks:
        velfile = os.path.join(vel_base_directory + network, 'ITRF_' + network + '_vels.txt');
        if os.path.isfile(velfile):
            name = np.loadtxt(velfile, skiprows=3, usecols=(0,), unpack=True, dtype={'names': ('name',),
                                                                                     'formats': ('U4',)})
            for station in name:
                download_usgs_ts_file(station[0], network, ts_base_directory)    
    return;


if __name__ == "__main__":
    # Go get all USGS GPS files from the website!
    # This gets velocity tables and time series in ITRF2008/NA
    # When you're done downloading, you must run the cache function
    
    download_usgs_velocities();
    download_usgs_offsets();
    download_usgs_time_series();
    cache_usgs_start_endtimes();
