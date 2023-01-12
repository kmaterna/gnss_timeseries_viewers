"""
Functions to generate a list of GNSS station names/locations within certain radius of a point or within a box.
num_years is minimum duration of time series
max_sigma is formal uncertainty on velocity in mm
"""

import numpy as np
import matplotlib.path as mpltPath
from Tectonic_Utils.geodesy import haversine
from Tectonic_Utils.read_write import read_kml
from . import vel_functions, gps_input_vel_pipeline


# DRIVER 1: STATIONS WITHIN RADIUS
def get_stations_within_radius(data_config_file, center, radius, coord_box=(-126, -110, 30.0, 49), network='pbo',
                               num_years=3.0, max_sigma=2, refframe="ITRF"):
    """
    :param data_config_file: string, filename
    :param center: list, [lon, lat]
    :param radius: float, km
    :param coord_box: tuple, (W, E, S, N), default is western US
    :param network: string, default 'pbo'
    :param num_years: float, default 3.0
    :param max_sigma: float, mm/yr, default 2.0
    :param refframe: string, default 'ITRF'
    :returns: list of strings, list of floats, list of floats, list of floats [stations, lon, lat, distances]
    """
    myVelfield = inputs_velfield(data_config_file, network, num_years, max_sigma, coord_box, refframe=refframe);
    close_stations, lon, lat, rad_distance = compute_circle(myVelfield, center, radius);
    return close_stations, lon, lat, rad_distance;


# DRIVER 2: STATIONS WITHIN BOX
def get_stations_within_box(data_config_file, coord_box, network='pbo', num_years=3.0, max_sigma=2, refframe="ITRF"):
    """
    :param data_config_file: string, filename
    :param coord_box: tuple, (W, E, S, N), required
    :param network: string, default 'pbo'
    :param num_years: float, default 3.0
    :param max_sigma: float, mm/yr, default 2.0
    :param refframe: string, default 'ITRF'
    :returns: list of strings, list of floats, list of floats [stations, lon, lat]
    """
    myVelfield = inputs_velfield(data_config_file, network, num_years, max_sigma, coord_box, refframe=refframe);
    within_stations, lon, lat = compute_box(myVelfield, coord_box);
    return within_stations, lon, lat;


# DRIVER 3: STATIONS WITHIN A USER-DEFINED KML BORDER
def get_stations_within_polygon(data_config_file, polygon_file, coord_box=(-126, -110, 30.0, 49), network='pbo',
                                num_years=3.0, max_sigma=2, refframe="ITRF"):
    """
    :param data_config_file: string, filename
    :param polygon_file: string, filename to simple KML file
    :param coord_box: tuple, (W, E, S, N), default western US
    :param network: string, default 'pbo'
    :param num_years: float, default 3.0
    :param max_sigma: float, mm/yr, default 2.0
    :param refframe: string, default 'ITRF'
    :returns: list of strings, list of floats, list of floats [stations, lon, lat]
    """
    [polygon_lon, polygon_lat] = read_kml.read_simple_kml(polygon_file);
    myVelfield = inputs_velfield(data_config_file, network, num_years, max_sigma, coord_box, refframe=refframe);
    within_stations, lon, lat = compute_within_polygon(myVelfield, polygon_lon, polygon_lat);
    return within_stations, lon, lat;


# ------------ INPUTS ------------------ #
def inputs_velfield(data_config_file, network, num_years, max_sigma, coord_box, refframe='ITRF'):
    # Purpose: generate input velocity field.
    print("Searching for GNSS stations in the %s network " % network);
    myVelfield = gps_input_vel_pipeline.import_velfield(data_config_file, network, refframe=refframe);
    myVelfield = vel_functions.clean_velfield(myVelfield, num_years=num_years, max_horiz_sigma=max_sigma,
                                              max_vert_sigma=max_sigma * 3, coord_box=coord_box);
    myVelfield = vel_functions.remove_duplicates(myVelfield);
    return myVelfield;


# ----------- CIRCLE FUNCTIONS ---------------- # 
def compute_circle(myVelfield, center, radius):
    """
    :param myVelfield: list of Station_Vels
    :param center: list [lon, lat]
    :param radius: float, km
    :returns: list of strings, list of floats, list of floats, list of floats [stations, lon, lat, distances]
    """
    close_stations, rad_distance = [], [];
    lon, lat = [], [];
    for station_vel in myVelfield:
        mydist = haversine.distance([center[1], center[0]], [station_vel.nlat, station_vel.elon])
        if mydist <= radius:
            rad_distance.append(mydist);
            lon.append(station_vel.elon);
            lat.append(station_vel.nlat);
            close_stations.append(station_vel.name);
    print("Returning %d stations within %.3f km of %.4f, %.4f" % (len(close_stations), radius, center[0], center[1]));
    return close_stations, lon, lat, rad_distance;


# ----------- BOX FUNCTIONS ---------------- #
def compute_box(myVelfield, coord_box):
    """
    :param myVelfield: list of Station_Vels
    :param coord_box: list [W, E, S, N]
    :returns: list of strings, list of floats, list of floats [stations, lon, lat]
    """
    close_stations, lon, lat = [], [], [];
    for station_vel in myVelfield:
        if coord_box[0] <= station_vel.elon <= coord_box[1]:
            if coord_box[2] <= station_vel.nlat <= coord_box[3]:
                close_stations.append(station_vel.name);
                lon.append(station_vel.elon);
                lat.append(station_vel.nlat);
    print("Returning %d stations within box" % (len(close_stations)));
    return close_stations, lon, lat;


# ----------- POLYGON FUNCTIONS ---------------- #
def compute_within_polygon(myVelfield, polygon_lon, polygon_lat):
    """
    :param myVelfield: velfield object
    :param polygon_lon: list of lon polygon vertices
    :param polygon_lat: list of lat polygon vertices
    :returns: list of strings, list of floats, list of floats [stations, lon, lat]
    """
    lon, lat = [], [];
    polygon = np.row_stack((np.array(polygon_lon), np.array(polygon_lat))).T;
    points = np.row_stack((np.array(myVelfield.elon), myVelfield.nlat)).T;
    path = mpltPath.Path(polygon);
    inside2 = path.contains_points(points);
    within_stations = [i for (i, v) in zip(myVelfield.name, inside2) if v];
    print("Returning %d stations within the given polygon" % (len(within_stations)));
    return within_stations, lon, lat;
