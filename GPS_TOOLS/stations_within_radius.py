# Python functions to take an input coordinate, and 
# determine which GNSS stations are within a certain radius of that point or a given box.

import numpy as np
import haversine
import read_kml
import gps_io_functions
import matplotlib.path as mpltPath


# Reference: Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);

# DRIVER 1: STATIONS WITHIN RADIUS
def get_stations_within_radius(center, radius, coord_box=(-126, -110, 30.0, 49), network='pbo'):
    # default coord box is western US, but it can be reduced
    [input_file, num_years, max_sigma] = general_config(network);
    myVelfield = inputs(input_file, num_years, max_sigma, coord_box, network);
    close_stations, lon, lat, rad_distance = compute_circle(myVelfield, center, radius);
    return close_stations, lon, lat, rad_distance;


# DRIVER 2: STATIONS WITHIN BOX
def get_stations_within_box(coord_box, network='pbo'):
    [input_file, num_years, max_sigma] = general_config(network);
    myVelfield = inputs(input_file, num_years, max_sigma, coord_box, network);
    within_stations, lon, lat = compute_box(myVelfield, coord_box);
    return within_stations, lon, lat;


# DRIVER 3: STATIONS WITHIN A USER-DEFINED KML BORDER
def get_stations_within_polygon(polygon_file, coord_box=(-126, -110, 30.0, 49), network='pbo'):
    [input_file, num_years, max_sigma] = general_config(network);
    [polygon_lon, polygon_lat] = read_kml.read_simple_kml(polygon_file);
    myVelfield = inputs(input_file, num_years, max_sigma, coord_box, network);
    within_stations, lon, lat = compute_within_polygon(myVelfield, polygon_lon, polygon_lat);
    return within_stations, lon, lat;


# ----------- CONFIGURE OPTIONS ---------- # 
def general_config(network):
    myParams = gps_io_functions.read_config_file();
    if network == 'pbo' or network == 'cwu' or network == 'nmt':
        input_file = myParams.pbo_velocities;
    elif network == 'unr':
        input_file = myParams.unr_velocities;
    else:
        print("ERROR: Network %s not recognized (try pbo/cwu/nmt/unr)" % network);
        input_file = None;
    num_years = 3.0;
    max_sigma = 2.0;
    return [input_file, num_years, max_sigma];


# ------------ INPUTS ------------------ #
def inputs(input_file, num_years, max_sigma, coord_box, network):
    # Purpose: generate input velocity field.
    if network == 'pbo' or network == 'cwu' or network == 'nmt':
        [myVelfield] = gps_io_functions.read_pbo_vel_file(input_file);
    else:  # network = unr
        [myVelfield] = gps_io_functions.read_unr_vel_file(input_file);
    [myVelfield] = gps_io_functions.clean_velfield(myVelfield, num_years, max_sigma, coord_box);
    [myVelfield] = gps_io_functions.remove_duplicates(myVelfield);
    return myVelfield;


# ----------- CIRCLE FUNCTIONS ---------------- # 
def compute_circle(myVelfield, center, radius):
    close_stations = [];
    rad_distance = [];
    lon = [];
    lat = [];
    for i in range(len(myVelfield.name)):
        mydist = haversine.distance([center[1], center[0]], [myVelfield.nlat[i], myVelfield.elon[i]])
        if mydist <= radius:
            rad_distance.append(mydist);
            lon.append(myVelfield.elon[i]);
            lat.append(myVelfield.nlat[i]);
            close_stations.append(myVelfield.name[i]);
    print("Returning %d stations that are within %f km of center %f, %f" % (len(close_stations), radius, center[0], center[1]))
    return close_stations, lon, lat, rad_distance;


# ----------- BOX FUNCTIONS ---------------- # 

def compute_box(myVelfield, coord_box):
    close_stations = [];
    lon = [];
    lat = [];
    for i in range(len(myVelfield.name)):
        if myVelfield.elon[i] >= coord_box[0] and myVelfield.elon[i] <= coord_box[1]:
            if myVelfield.nlat[i] >= coord_box[2] and myVelfield.elon[i] <= coord_box[3]:
                close_stations.append(myVelfield.name[i]);
                lon.append(myVelfield.elon[i]);
                lat.append(myVelfield.nlat[i]);
    print("Returning %d stations that are within box" % (len(close_stations)));
    return close_stations, lon, lat;


# ----------- POLYGON FUNCTIONS ---------------- # 

def compute_within_polygon(myVelfield, polygon_lon, polygon_lat):
    # If within the polygon:
    lon = [];
    lat = [];  # not currently coded
    polygon = np.row_stack((np.array(polygon_lon), np.array(polygon_lat))).T;
    points = np.row_stack((np.array(myVelfield.elon), myVelfield.nlat)).T;
    path = mpltPath.Path(polygon);
    inside2 = path.contains_points(points);
    within_stations = [i for (i, v) in zip(myVelfield.name, inside2) if v];
    print("Returning %d stations that are within the given polygon" % (len(within_stations)));
    return within_stations, lon, lat;


if __name__ == "__main__":
    center = [-124, 40];
    radius = 120;
    close_stations, lons, lats, rad_distance = get_stations_within_radius(center, radius);
