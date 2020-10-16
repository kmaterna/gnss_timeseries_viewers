import requests
import pandas as pd
import datetime as dt

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


def download_usgs_velocity_tables(network, base_directory, outfile):
    url = base_url + network + '/velocities';
    df_list = pd.read_html(requests.get(url).content)
    print("Writing %d tables on website %s " % (len(df_list), url) );
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

        # Filtered
        filt_outfile = base_directory + "FILT_" + outfile;
        one_outfile = open(filt_outfile, 'w');
        one_outfile.write("# USGS Velocity File, Filtered\n")
        one_outfile.write("# Downloaded from %s on %s\n" % (url, dt.datetime.now()));
        one_outfile.write("# ");
        one_outfile.close();
        df_list[2].to_csv(filt_outfile, index=False, sep=' ', mode='a');
        print("Writing %s " % filt_outfile);
    return;


if __name__ == "__main__":
    base_directory = "../../GPS_POS_DATA/USGS_Data/Velocities/"
    download_usgs_velocity_tables("NCalifornia_SGPS", base_directory, "NCalifornia_SGPS_vels.txt");
    download_usgs_velocity_tables("Pacific_Northwest", base_directory, "Pacific_Northwest_vels.txt");
    download_usgs_velocity_tables("SFBayArea", base_directory, "SFBayArea_vels.txt");
    download_usgs_velocity_tables("SFBayArea_SGPS", base_directory, "SFBayArea_SGPS_vels.txt");

