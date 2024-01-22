from gnss_timeseries_viewers.gps_tools import load_gnss, offsets, single_station_tsplot
import datetime as dt
import unittest, os

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))  # 4 dirs up
config_file = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt');

class Tests(unittest.TestCase):
    def test_read_unr_ts(self):
        station = "P499"
        database = load_gnss.create_station_repo(config_file, proc_center='unr', refframe='NA');
        [myData, offset_obj, eq_obj] = database.load_station(station);

        # Remove outliers
        myData = myData.remove_outliers(outliers_def=20)

        # Impose time limits
        myData = myData.impose_time_limits(dt.datetime.strptime("2000-01-01", "%Y-%m-%d"),
                                           dt.datetime.strptime("2025-01-01", "%Y-%m-%d"))
        # Remove offsets and antenna changes
        nooffsets = offsets.remove_offsets(myData, offset_obj)

        # Remove earthquakes
        nooffsets = offsets.remove_offsets(nooffsets, eq_obj)

        single_station_tsplot.single_ts_plot(myData, detrended=nooffsets, savename=station+'_view_ts.png',
                                             title='GNSS time series from station '+station,
                                             detrended_label='no offsets')


if __name__ == "__main__":
    unittest.main()
