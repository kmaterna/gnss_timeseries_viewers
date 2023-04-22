"""
Functions to operate on Station-Vel objects, or velocity-field objects (lists of station-vels)
Future work: Combining two velocity fields in outer combination (inclusive) or inner combination (common stations)
"""
from . import utilities
import numpy as np
import matplotlib.path as mpltPath
from Tectonic_Utils.geodesy import haversine
import Tectonic_Utils.geodesy.xyz2llh as geo_conv
from Tectonic_Utils.geodesy import fault_vector_functions


class Station_Vel:
    def __init__(self, name, nlat, elon, n, e, u, sn, se, su, first_epoch, last_epoch,
                 refframe=None, proccenter=None, subnetwork=None, survey=None, meas_type='gnss'):
        self.name = name;
        self.nlat = nlat;
        self.elon = utilities.check_lon_sanity(elon);  # -180 < lon < 180, used for velfields
        self.n = n;  # in mm/yr
        self.e = e;  # in mm/yr
        self.u = u;  # in mm/yr
        self.sn = sn;
        self.se = se;
        self.su = su;
        self.first_epoch = first_epoch;
        self.last_epoch = last_epoch;
        self.refframe = refframe;
        self.proccenter = proccenter;
        self.subnetwork = subnetwork;
        self.survey = survey;
        self.meas_type = meas_type;

    # ---------- PREDICATES ------------- #

    def spans_shorter_than(self, num_years):
        """
        :param num_years: float, number of years
        :return: bool
        """
        deltatime = self.last_epoch - self.first_epoch;
        if deltatime.days <= num_years * 365.24:
            return 1;   # too short time interval
        else:
            return 0;

    def spans_longer_than(self, num_years):
        """
        :param num_years: float, number of years
        :return: bool
        """
        deltatime = self.last_epoch - self.first_epoch;
        if deltatime.days >= num_years * 365.24:
            return 1;
        else:
            return 0;   # too short time interval

    def within_bbox(self, bbox):
        """
        :param bbox: [W, E, S, N] list of 4 floats
        :return: bool
        """
        if bbox[0] <= self.elon <= bbox[1] and bbox[2] <= self.nlat <= bbox[3]:
            return 1;
        else:
            return 0;

    def survives_blacklist(self, blacklist):
        """
        :param blacklist: list of strings
        :return: bool
        """
        if self.name in blacklist:
            return 0;
        else:
            return 1;

    def get_azimuth_degrees(self):
        """
        :return: azimuth of velocity or offset, in degrees clockwise from north
        """
        if fault_vector_functions.get_vector_magnitude([self.e, self.n]) <= 0.000001:
            azimuth = 0;
        else:
            azimuth = fault_vector_functions.get_strike(self.e, self.n);
        return azimuth;


class Station_Vel_XYZ:
    # Station_Vel_XYZ are in m/yr (under consideration), with XYZ as ECEF position in meters. Used for velocities.
    def __init__(self, name, x_pos, y_pos, z_pos, x_rate, y_rate, z_rate, x_sigma, y_sigma, z_sigma,
                 first_epoch, last_epoch):
        self.name = name;
        self.x_pos = x_pos;
        self.y_pos = y_pos;
        self.z_pos = z_pos;
        self.x_rate = x_rate;
        self.y_rate = y_rate;
        self.z_rate = z_rate;
        self.x_sigma = x_sigma;
        self.y_sigma = y_sigma;
        self.z_sigma = z_sigma;
        self.first_epoch = first_epoch;
        self.last_epoch = last_epoch;


# ---------- BASIC FUNCTIONS ON VELFIELDS ------------- #

def basic_clean_stations(velfield, coord_box=(-180, 180, -90, 90), num_years=3.0, max_sigma=2):
    velfield = clean_velfield(velfield, num_years, max_sigma, max_sigma * 3, coord_box);
    velfield = remove_duplicates(velfield);
    return velfield;


def clean_velfield(velfield, num_years=0.0, max_horiz_sigma=1000, max_vert_sigma=1000, coord_box=(-180, 180, -90, 90),
                   verbose=False):
    """
    Filter into cleaner GPS velocities by time series length, formal uncertainty, or geographic range.
    Default arguments are meant to be global.
    Verbose means print why you're excluding stations.

    :param velfield: a list of station_vel objects
    :param num_years: number of years, float
    :param max_horiz_sigma: formal uncertainty for east or north, mm/yr
    :param max_vert_sigma: formal uncertainty for vertical, mm/yr
    :param coord_box: list or tuple geographic range, [w, e, s, n]
    :param verbose: bool
    """
    cleaned_velfield = [];
    if verbose:
        print("Cleaning: Starting with " + str(len(velfield)) + " stations");
    for station in velfield:
        if station.sn > max_horiz_sigma:  # too high sigma, please exclude
            print('Excluding ' + station.name + ' for large north sigma of ' + str(
                station.sn) + ' mm/yr') if verbose else 0;
            continue;
        if station.se > max_horiz_sigma:
            print('Excluding ' + station.name + ' for large east sigma of ' + str(
                station.se) + ' mm/yr') if verbose else 0;
            continue;
        if station.su > max_vert_sigma:
            print('Excluding ' + station.name + ' for large vertical sigma of' + str(
                station.su) + ' mm/yr') if verbose else 0;
            continue;
        if station.spans_shorter_than(num_years):  # too short time interval, please exclude
            print('Excluding ' + station.name + 'for time range too short') if verbose else 0;
            continue;
        if station.within_bbox(coord_box):  # The station is within the box of interest.
            cleaned_velfield.append(station);
        else:
            print("excluding for outside box of interest: ", station.elon, station.nlat) if verbose else 0;
    if verbose:
        print("Cleaning: After applying selection criteria, we have " + str(len(cleaned_velfield)) + " stations\n");
    return cleaned_velfield;


def remove_duplicates(velfield, verbose=False):
    """
    Right now assuming all entries at the same lat/lon are same station
    """
    cleaned_velfield = []
    if verbose:
        print("Removing duplicates: Starting with %d stations" % len(velfield));
    for vel in velfield:
        is_duplicate = 0;
        for comp_station_vel in cleaned_velfield:
            if abs(comp_station_vel.elon - vel.elon) < 0.0005 and abs(comp_station_vel.nlat - vel.nlat) < 0.0005:
                # we found a duplicate measurement.
                is_duplicate = 1;
        if is_duplicate == 0:
            cleaned_velfield.append(vel);
    if verbose:
        print("Removing duplicates: Ending with %d stations" % len(cleaned_velfield));
    return cleaned_velfield;


def disp_points_to_station_vels(obs_disp_points):
    """
    Convert from disp_point_object (which might contain velocities) to station_vel object.
    This is the opposite of Elastic_stresses_py.PyCoulomb.disp_points_object.utilities.station_vel_object_to_disp_points
    """
    station_vel_list = [];
    for item in obs_disp_points:
        new_obj = Station_Vel(nlat=item.lat, elon=item.lon, name=item.name, n=item.dN_obs * 1000,
                              e=item.dE_obs * 1000, u=item.dU_obs * 1000, sn=item.Sn_obs * 1000,
                              se=item.Se_obs * 1000, su=item.Su_obs * 1000, meas_type=item.meas_type,
                              first_epoch=item.starttime, last_epoch=item.endtime, refframe=item.refframe);
        station_vel_list.append(new_obj);
    return station_vel_list;


def remove_blacklist_paired_data(velfield, blacklist, matching_data=None, verbose=False):   # ZIP VERSION
    """
    :param velfield: list of StationVels
    :param blacklist: list of strings
    :param verbose: boolean, false
    :param matching_data: a vector of things that go with each station, like epicentral distance
    :return: list of StationVels, list of matching data
    """
    cleaned_velfield, cleaned_distances = [], [];
    if verbose:
        print("Removing blacklist: Starting with %d stations" % len(velfield));
    for i in range(len(velfield)):
        if velfield[i].survives_blacklist(blacklist):
            cleaned_velfield.append(velfield[i]);
            cleaned_distances.append(matching_data[i]);
    if verbose:
        print("Removing blacklist: Ending with %d stations" % len(cleaned_velfield));
    return cleaned_velfield, cleaned_distances;


def remove_blacklist_vels(velfield, blacklist, verbose=False):   # NON-ZIP VERSION
    """
    :param velfield: list of StationVels
    :param blacklist: list of strings
    :param verbose: boolean, false
    :return: list of StationVels
    """
    cleaned_velfield = [];
    if verbose:
        print("Removing blacklist: Starting with %d stations" % len(velfield));
    for i in range(len(velfield)):
        if velfield[i].survives_blacklist(blacklist):
            cleaned_velfield.append(velfield[i]);
    if verbose:
        print("Removing blacklist: Ending with %d stations" % len(cleaned_velfield));
    return cleaned_velfield;


def filter_to_circle(myVelfield, center, radius):
    """
    :param myVelfield: list of Station_Vels
    :param center: list [lon, lat]
    :param radius: float, km
    :returns: list of Station_Vels, list of floats [stations, distances]
    """
    close_stations, rad_distance = [], [];
    for station_vel in myVelfield:
        mydist = haversine.distance([center[1], center[0]], [station_vel.nlat, station_vel.elon])
        if mydist <= radius:
            rad_distance.append(mydist);
            close_stations.append(station_vel);
    print("Returning %d stations within %.3f km of %.4f, %.4f" % (len(close_stations), radius, center[0], center[1]));
    return close_stations, rad_distance;


def filter_to_bounding_box(myVelfield, coord_box):
    """
    :param myVelfield: list of Station_Vels
    :param coord_box: list [W, E, S, N]
    :returns: list of station_vels
    """
    close_stations = [];
    for station_vel in myVelfield:
        if station_vel.within_bbox(coord_box):
            close_stations.append(station_vel);
    print("Returning %d stations within box" % (len(close_stations)));
    return close_stations;


def filter_to_within_polygon(myVelfield, polygon_lon, polygon_lat):
    """
    :param myVelfield: velfield object
    :param polygon_lon: list of lon polygon vertices
    :param polygon_lat: list of lat polygon vertices
    :returns: list of station_vels
    """
    polygon = np.row_stack((np.array(polygon_lon), np.array(polygon_lat))).T;
    points = np.row_stack((np.array(myVelfield.elon), myVelfield.nlat)).T;
    path = mpltPath.Path(polygon);
    inside2 = path.contains_points(points);
    within_stations = [i for (i, v) in zip(myVelfield, inside2) if v];
    print("Returning %d stations within the given polygon" % (len(within_stations)));
    return within_stations;


def velocity_misfit_function(Velfield1, Velfield2):
    """
    Compares two velocity fields with identical list of stations by taking RMS of velocity differences.
    """
    misfit = 0;
    for i in range(len(Velfield1)):
        misfit += np.square(Velfield2[i].e - Velfield1[i].e);
        misfit += np.square(Velfield2[i].n - Velfield1[i].n);
        misfit += np.square(Velfield2[i].u - Velfield1[i].u);
    misfit /= 3 * len(Velfield1);
    misfit = np.sqrt(misfit);
    return misfit;


def convert_enu_velfield_to_xyz(velfield):
    """
    Convert Station_Vel objects into Station_Vel_XYZ objects, with position/velocity/unc in ECEF-XYZ coordinates.
    All units are in meters or meters/year
    Right now, there's no way to bring elevation information into this calculation.
    """
    xyz_objects = [];
    for item in velfield:
        llh_origin = np.array([[item.elon, item.nlat, 0], ]);  # position vector lon, lat
        llh_origin_simple = np.array([item.elon, item.nlat, 0]);  # simpler numpy array
        enu_stds = 0.001 * np.array([item.se, item.sn, item.su]);
        xyz_pos = geo_conv.llh2xyz(llh_origin);
        enu_vel = 0.001 * np.array([[item.e, item.n, item.u], ]);  # velocity in m
        ecov = np.diag(enu_stds);  # a 3x3 np array matrix with covariances on the diagonal
        xyz_vel, xyz_cov = geo_conv.enu2xyz(enu_vel, llh_origin_simple, ecov);
        newobj = Station_Vel_XYZ(name=item.name, x_pos=xyz_pos[0][0], y_pos=xyz_pos[0][1],
                                 z_pos=xyz_pos[0][2], x_rate=xyz_vel[0][0], y_rate=xyz_vel[0][1],
                                 z_rate=xyz_vel[0][2], x_sigma=xyz_cov[0][0], y_sigma=xyz_cov[1][1],
                                 z_sigma=xyz_cov[2][2], first_epoch=item.first_epoch,
                                 last_epoch=item.last_epoch);
        xyz_objects.append(newobj);
    return xyz_objects;


def get_bounding_box(velfield, border=0.1):
    """
    Get a bounding box for a list of station_vels, padded by a border region.

    :param velfield: list of StationVels
    :param border: float, padding for bounding box, in degrees
    :returns: list of [W, E, S, N]
    """
    lons = [x.elon for x in velfield]
    lats = [x.nlat for x in velfield]
    bbox = [np.min(lons) - border, np.max(lons) + border, np.min(lats) - border, np.max(lats) + border];
    return bbox;
