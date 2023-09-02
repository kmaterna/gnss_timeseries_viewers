
import unittest
import GNSS_TimeSeries_Viewers.gps_tools.load_gnss as load_gnss

config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"

class Tests(unittest.TestCase):

    def test_unr_database(self):
        database = load_gnss.create_station_repo(config_file, proc_center='unr', refframe='NA');
        database.import_full_velfield();
        [_myData, _, _] = database.load_station("P496");
        database = load_gnss.create_station_repo(config_file, proc_center='unr', refframe='ITRF');
        # database.import_full_velfield();  # TODO
        [_myData, _, _] = database.load_station("P496");

    def test_cwu_database(self):
        database = load_gnss.create_station_repo(config_file, proc_center='cwu', refframe='NA');
        database.import_full_velfield();
        [_myData, _, _] = database.load_station("P496");
        database = load_gnss.create_station_repo(config_file, proc_center='cwu', refframe='ITRF');
        database.import_full_velfield();
        [myData, _, _] = database.load_station("P496");
        myData = myData.remove_nans();  # testing a computational function
        _slopes = myData.get_slope();

    def test_pbo_database(self):
        database = load_gnss.create_station_repo(config_file, proc_center='pbo', refframe='NA');
        database.import_full_velfield();
        [_myData, _, _] = database.load_station("P496");
        database = load_gnss.create_station_repo(config_file, proc_center='pbo', refframe='ITRF');
        database.import_full_velfield();
        [_myData, _, _] = database.load_station("P496");

    def test_usgs_database(self):
        database = load_gnss.create_station_repo(config_file, proc_center='usgs', refframe='NA',
                                                 subnetwork='SFBayArea');
        database.import_full_velfield();
        [_myData, _, _] = database.load_station("P236");
        database = load_gnss.create_station_repo(config_file, proc_center='usgs', refframe='ITRF',
                                                 subnetwork='SFBayArea');
        [_myData, _, _] = database.load_station("P236");
        database.import_full_velfield();


if __name__ == "__main__":
    unittest.main()
