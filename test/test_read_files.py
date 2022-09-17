
import unittest
import GNSS_TimeSeries_Viewers.gps_tools as gps_tools

config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"

class Tests(unittest.TestCase):

    def test_reading_ts(self):
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "unr", config_file, "NA");
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "unr", config_file, "ITRF");
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "pbo", config_file, "NA");
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "pbo", config_file, "ITRF");
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "cwu", config_file, "NA");
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "cwu", config_file, "ITRF");
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "usgs", config_file, "NA");
        [_myData, _, _] = gps_tools.gps_input_pipeline.get_station_data("P496", "usgs", config_file, "ITRF");

    def test_reading_velocities(self):
        _v = gps_tools.gps_input_vel_pipeline.import_velfield(config_file, network='pbo', refframe='ITRF');
        _v = gps_tools.gps_input_vel_pipeline.import_velfield(config_file, network='unr', refframe='NA');
        _v = gps_tools.gps_input_vel_pipeline.import_velfield(config_file, network='cwu', refframe='NA');
        _v = gps_tools.gps_input_vel_pipeline.import_velfield(config_file, network='usgs', refframe='NA',
                                                              sub_network='SFBayArea');
        _v = gps_tools.gps_input_vel_pipeline.import_velfield(config_file, network='usgs', refframe='ITRF',
                                                              sub_network='SFBayArea');


if __name__ == "__main__":
    unittest.main()
