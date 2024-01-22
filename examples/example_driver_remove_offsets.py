#!/usr/bin/env python

"""
Give an example workflow of removing offsets from a time series, both through default functions and
through custom techniques for the Ridgecrest earthquake scenario.
The example demonstrates a few other time series functions as well, such as outliers and time limits.

Kathryn Materna
1/22/2024
"""

from gnss_timeseries_viewers.gps_tools import load_gnss, offsets, single_station_tsplot
import datetime as dt
import os

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # 4 dirs up
data_config_file = os.path.join(base_dir, 'GEOPHYS_DATA', 'GPS_POS_DATA', 'config.txt')
ridgecrest_window = [dt.datetime.strptime("2019-07-04", "%Y-%m-%d"), dt.datetime.strptime("2019-07-07", "%Y-%m-%d")]


def load_station_remove_offsets():
    station = "P595"
    database = load_gnss.create_station_repo(data_config_file, proc_center='unr', refframe='NA');  # point to configfile
    [myData, offset_obj, eq_obj] = database.load_station(station);  # gives Data, antenna changes, and eq steps

    # Example function: Impose time limits
    myData = myData.impose_time_limits(dt.datetime.strptime("2000-01-01", "%Y-%m-%d"),
                                       dt.datetime.strptime("2025-01-01", "%Y-%m-%d"))

    # Remove outliers > 20 mm for clarity
    myData = myData.remove_outliers(outliers_def=20)

    # Remove offsets due to antenna changes
    nooffsets = offsets.remove_offsets(myData, offset_obj)

    # Remove all earthquake offsets, the normal way for most events and the special way for Ridgecrest
    for single_earthquake_offset in eq_obj:

        if ridgecrest_window[0] < single_earthquake_offset.evdt < ridgecrest_window[1]:  # the special way
            print("Special window for Ridgecrest, which had 2 M>6.4 earthquakes in 2 days")
            ridgecrest_coseismic = offsets.solve_for_offsets(nooffsets, [ridgecrest_window], num_days=10)
            nooffsets = offsets.remove_offsets(nooffsets, ridgecrest_coseismic)
            # Warning: Some pathological values may appear on 7/5/19 and 7/6/19 since they're within the solving window
            # You might want to remove them with the following function:
            nooffsets = nooffsets.remove_specific_date(dt.datetime.strptime("2019-07-05", "%Y-%m-%d"))
            nooffsets = nooffsets.remove_specific_date(dt.datetime.strptime("2019-07-06", "%Y-%m-%d"))

        else:  # remove offsets the normal way
            nooffsets = offsets.remove_offsets(nooffsets, [single_earthquake_offset])  # normal way, use list of offsets

    single_station_tsplot.single_ts_plot(myData, detrended=nooffsets, savename=station+'_view_ts.png',
                                         title='GNSS time series from station '+station, outdir='example_pngs/',
                                         detrended_label='no offsets')  # simple before-and-after plot


if __name__ == "__main__":
    load_station_remove_offsets()
