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
vel_base_directory = "../../GPS_POS_DATA/USGS_Data/Velocities/"
ts_base_directory = "../../GPS_POS_DATA/USGS_DATA/Time_Series/"
usgs_cache_file = "../../GPS_POS_DATA/USGS_DATA/usgs_station_cache.txt"


def download_usgs_velocity_tables(network, base_directory, outfile):
    # Go online and grab USGS velocity tables, in NA, ITRF2008, and Filtered
    # One table for each network, as organized in the USGS database.
    url = base_url + network + '/velocities';
    df_list = pd.read_html(requests.get(url).content)
    print("Reading %d tables on website %s " % (len(df_list), url) );
    if len(df_list) == 3:
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
        print("Skipping filtered velocity field for convenience.")
    else:
        print("Skipping network %s because of unexpected number of reference frames." % network)
    return;


def download_usgs_time_series(station, network, ts_base_directory):
    # Download some time series files from USGS in NA and ITRF2008
    # They will be organized by network, as they are in the USGS database
    directory_name = ts_base_directory + network + '/';
    subprocess.call(['mkdir', '-p', directory_name], shell=False);
    url_name = base_url + 'data/networks/' + network + '/' + station.lower() + '/nafixed/' + station.lower() + '.rneu'
    subprocess.call(["wget", url_name, '-O', directory_name+station.lower()+'_NAfixed.rneu'], shell=False);
    url_name = base_url + 'data/networks/' + network + '/' + station.lower() + '/itrf2008/' + station.lower() + '.rneu'
    subprocess.call(["wget", url_name, '-O', directory_name+station.lower()+'_ITRF2008.rneu'], shell=False);
    return;


def download_usgs_vels_and_ts():
    for network in networks:
        download_usgs_velocity_tables(network, vel_base_directory, outfile=network+"_vels.txt");
        velfile = vel_base_directory+'ITRF_'+network+'_vels.txt';
        if os.path.isfile(velfile):
            name = np.loadtxt(velfile, skiprows=3, usecols=(0,), unpack=True, dtype={'names': ('name',), 'formats': ('U4',)})
            for station in name:
                download_usgs_time_series(station[0], network, ts_base_directory)
    return;


def usgs_network_directory_from_velfile(infile):
    # Network parsing and plumbing
    # We assume a parallel directory above the Velocity directory with a bunch of TS files
    # From which we can derive startdate and enddate
    usgs_network = infile.split('/')[-1][0:-9];
    print(usgs_network);
    if usgs_network[0:4] == 'NAM_':
        usgs_network = usgs_network[4:];
    if usgs_network[0:5] == 'ITRF_':
        usgs_network = usgs_network[5:];
    usgs_directory = '';
    for i in range(len(infile.split('/')) - 2):
        usgs_directory = usgs_directory + infile.split('/')[i];
        usgs_directory = usgs_directory + '/';
    network_ts_directory = usgs_directory + 'Time_Series/' + usgs_network + '/'
    total_vel_directory = usgs_directory + 'Velocities/'
    return network_ts_directory, total_vel_directory;


def cache_usgs_vels_and_ts():
    print("Writing network cache %s " % (usgs_cache_file) );
    ofile = open(usgs_cache_file, 'w');
    for network in networks:
        velfile = vel_base_directory + 'ITRF_' + network + '_vels.txt';
        network_ts_directory, _ = usgs_network_directory_from_velfile(velfile);

        [names, lon, lat] = np.loadtxt(velfile, skiprows=3, usecols=(0, 1, 2), unpack=True, dtype={'names': (
                'name', 'lon', 'lat'), 'formats': ('U4', np.float, np.float)});
        # Populating the first_epoch and last_epoch with information from the associated time series directory.
        for i in range(len(names)):
            ts_filename = network_ts_directory + '/' + names[i] + '_NAfixed.rneu';
            dates = np.loadtxt(ts_filename, unpack=True, usecols=(0,),
                               dtype={'names': ('dtstrs',), 'formats': ('U8',)}, ndmin=1);
            first_epoch=dates[0][0];
            last_epoch=dates[-1][0];
            ofile.write("%s %.5f %.5f %s %s %s\n" % (names[i], lon[i], lat[i], first_epoch, last_epoch, network) );

    ofile.close();
    return;


if __name__ == "__main__":
    # Go get all the USGS GPS files from the internet!
    # This gets velocity tables from website, AND time series in ITRF2008/NA
    # When you're done downloading, you must run the cache function
    download_usgs_vels_and_ts();
    cache_usgs_vels_and_ts();
