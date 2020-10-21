import requests
import pandas as pd
import datetime as dt
import os
import subprocess
import gps_io_functions

# GLOBAL VARIABLES
base_url = 'https://earthquake.usgs.gov/monitoring/gps/'
networks = ['Alaska',
            'Alaska_SGPS',
            'BasinAndRange',
            'BasinAndRange_SGPS',
            'CentralCalifornia',
            'CentralCalifornia_ITRF2014',
            'CentralCalifornia_SGPS',
            'CentralUS',
            'DiabloCanyon_SPGPS',
            'EasternORWA_SGPS',
            'ECSZ_SGPS',
            'ElCentro',
            'Helens',
            'JuanDeFuca',
            'LKY2',
            'LongValley',
            'Mammoth_SGPS',
            'MtBaker',
            'NCalifornia_SGPS',
            'NorthEasternCal_SGPS',
            'Northridge',
            'Pacific',
            'Pacific_Northwest',
            'SFBayArea',
            'SFBayArea_SGPS',
            'Sisters',
            'Southern_California',
            'WindKetchFlat_SGPS',
            'YellowstoneContin',
            'Yellowstone_SPGPS'];
vel_base_directory = "../../GPS_POS_DATA/USGS_Data/Velocities/"
ts_base_directory = "../../GPS_POS_DATA/USGS_DATA/Time_Series/"


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
    subprocess.call(["wget", url_name,'-O',directory_name+station.lower()+'_NAfixed.rneu'],shell=False);
    url_name = base_url + 'data/networks/' + network + '/' + station.lower() + '/itrf2008/' + station.lower() + '.rneu'
    subprocess.call(["wget", url_name,'-O',directory_name+station.lower()+'_ITRF2008.rneu'],shell=False);
    return;


if __name__ == "__main__":
    # Go get all the USGS GPS files!
    # This gets most of them (anything with 3 tables: NA, ITRF, Filtered)
    # Velocity tables AND time series in ITRF2008/NA
    for network in networks:
        download_usgs_velocity_tables(network, vel_base_directory, outfile=network+"_vels.txt");
        velfile = vel_base_directory+'ITRF_'+network+'_vels.txt';
        ts_directory = ts_base_directory+network+'/';
        if network in ['Southern_California', 'WindKetchFlat_SGPS', 'YellowstoneContin','Yellowstone_SPGPS']:
            if os.path.isfile(velfile):
                [myVelfield] = gps_io_functions.read_usgs_velfile(velfile);
                for station in myVelfield.name:
                    download_usgs_time_series(station, network, ts_base_directory)
