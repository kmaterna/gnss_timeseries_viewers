"""
Functions to operate on velocity-field objects (lists of station-vels)
Future work: Combining two velocity fields in outer combination (inclusive) or inner combination (common stations)
"""

import numpy as np
import scipy.optimize
from . import gps_io_functions
import Tectonic_Utils.geodesy.xyz2llh as geo_conv


def clean_velfield(velfield, num_years=0, max_horiz_sigma=1000, max_vert_sigma=1000, coord_box=(-180, 180, -90, 90),
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
            print('Excluding '+station.name+' for large north sigma of '+str(station.sn)+' mm/yr') if verbose else 0;
            continue;
        if station.se > max_horiz_sigma:
            print('Excluding '+station.name+' for large east sigma of '+str(station.se)+' mm/yr') if verbose else 0;
            continue;
        if station.su > max_vert_sigma:
            print('Excluding '+station.name+' for large vertical sigma of'+str(station.su)+' mm/yr') if verbose else 0;
            continue;
        deltatime = station.last_epoch - station.first_epoch;
        if deltatime.days <= num_years * 365.24:  # too short time interval, please exclude
            print('Excluding ' + station.name + 'for time range too short') if verbose else 0;
            continue;
        if coord_box[0] < station.elon < coord_box[1] and coord_box[2] < station.nlat < coord_box[3]:
            # The station is within the box of interest.
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
        new_obj = gps_io_functions.Station_Vel(nlat=item.lat, elon=item.lon, name=item.name, n=item.dN_obs*1000,
                                               e=item.dE_obs*1000, u=item.dU_obs*1000, sn=item.Sn_obs*1000,
                                               se=item.Se_obs*1000, su=item.Su_obs*1000, meas_type=item.meas_type,
                                               first_epoch=item.starttime, last_epoch=item.endtime,
                                               refframe=item.refframe);
        station_vel_list.append(new_obj);
    return station_vel_list;


def remove_blacklist_vels(velfield, blacklist, verbose=False):
    cleaned_velfield = []
    if verbose:
        print("Removing blacklist: Starting with %d stations" % len(velfield));
    for vel in velfield:
        if vel.name not in blacklist:
            cleaned_velfield.append(vel);
    if verbose:
        print("Removing blacklist: Ending with %d stations" % len(cleaned_velfield));
    return cleaned_velfield;


def velocity_misfit_function(Velfield1, Velfield2):
    """
    Compares two velocity fields with identical list of stations by taking RMS of velocity differences.
    """
    misfit = 0;
    for i in range(len(Velfield1)):
        misfit += np.square(Velfield2[i].e - Velfield1[i].e);
        misfit += np.square(Velfield2[i].n - Velfield1[i].n);
        misfit += np.square(Velfield2[i].u - Velfield1[i].u);
    misfit /= 3*len(Velfield1);
    misfit = np.sqrt(misfit);
    return misfit;


def Apply_Helmert_Transformation(xyz_velfieldA, Helmert_params):
    """
    Apply a seven-parameter Helmert transformation with a given vector of seven parameters
    https://en.wikipedia.org/wiki/Helmert_transformation
    Only transforms positions (so we do something clever beforehand for velocities)

    :param xyz_velfieldA: Station_Vel_XZY objects in un-primed reference frame, will be transformed by Helmert trans.
    :type xyz_velfieldA: list
    :param Helmert_params: seven params, [Cx, Cy, Cz, s, Rx, Ry, Rz]
    :type Helmert_params: list
    :returns: xyz_velfieldB, list of Station_Vel_XZY objects in primed reference frame
    """
    xyz_velfieldB = [];
    [Cx, Cy, Cz, s, Rx, Ry, Rz] = [*Helmert_params];
    s_mult = 1 + s*1e-6;
    H_Matrix = np.array([[1, -Rz, Ry],
                         [Rz, 1, -Rx],
                         [-Ry, Rx, 1]]);
    scaled_H_Matrix = np.multiply(H_Matrix, s_mult);
    translation_vector = np.array([[Cx, Cy, Cz]]);
    for item in xyz_velfieldA:
        # We store raw coordinates in x_pos, and velocity-considered coordinates in x_rate
        position_in_A = np.array([[item.x_rate, item.y_rate, item.z_rate]]);   # We use rate for transformed position
        position_in_B = np.dot(scaled_H_Matrix, position_in_A.T);   # Performing rotation/scaling transformation
        position_in_B = np.add(position_in_B, translation_vector.T);  # adding translation transformation
        position_in_B = position_in_B.T;
        xyz_stds = np.array([item.x_sigma, item.y_sigma, item.z_sigma]);
        ecov = np.diag(xyz_stds);  # a 3x3 np array matrix with covariances on the diagonal
        sigma_B = np.dot(np.dot(scaled_H_Matrix, ecov), scaled_H_Matrix.T);   # Transform sigmas
        sigma_B = np.multiply(sigma_B, np.square(s_mult));
        newobj = gps_io_functions.Station_Vel_XYZ(name=item.name, x_pos=item.x_pos, y_pos=item.y_pos, z_pos=item.z_pos,
                                                  x_rate=position_in_B[0][0],  y_rate=position_in_B[0][1],
                                                  z_rate=position_in_B[0][2],
                                                  x_sigma=sigma_B[0][0], y_sigma=sigma_B[1][1], z_sigma=sigma_B[2][2],
                                                  first_epoch=item.first_epoch, last_epoch=item.last_epoch);
        xyz_velfieldB.append(newobj);
    return xyz_velfieldB;


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
        enu_vel = 0.001 * np.array([[item.e, item.n, item.u], ]);   # velocity in m
        ecov = np.diag(enu_stds);  # a 3x3 np array matrix with covariances on the diagonal
        xyz_vel, xyz_cov = geo_conv.enu2xyz(enu_vel, llh_origin_simple, ecov);
        newobj = gps_io_functions.Station_Vel_XYZ(name=item.name, x_pos=xyz_pos[0][0], y_pos=xyz_pos[0][1],
                                                  z_pos=xyz_pos[0][2], x_rate=xyz_vel[0][0], y_rate=xyz_vel[0][1],
                                                  z_rate=xyz_vel[0][2], x_sigma=xyz_cov[0][0], y_sigma=xyz_cov[1][1],
                                                  z_sigma=xyz_cov[2][2], first_epoch=item.first_epoch,
                                                  last_epoch=item.last_epoch);
        xyz_objects.append(newobj);
    return xyz_objects;


def prepare_velocities_for_helmert_trans(velfield):
    """
    Take regular Station_Vel objects (ENU).
    Then, we are using the "Station_Vel_XYZ" object in a clever way.
    x_pos, y_pos, z_pos are related to ONLY the lon/lat position
    x_rate, y_rate, z_rate are the position when a few years of velocity has been added
    x_sigma, y_sigma, z_sigma remain unchanged.
    This is going to be used for helmert transformation of positions into a new reference frame.
    """
    velfield_xyz = convert_enu_velfield_to_xyz(velfield);  # convert to xyz
    special_pos_objects = [];
    multiplier = 100;  # number of years for velocities, to avoid floating point and rounding errors
    for item in velfield_xyz:
        x_special_pos = item.x_pos + multiplier * item.x_rate;
        y_special_pos = item.y_pos + multiplier * item.y_rate;
        z_special_pos = item.z_pos + multiplier * item.z_rate;
        newobj = gps_io_functions.Station_Vel_XYZ(name=item.name, x_pos=item.x_pos, y_pos=item.y_pos, z_pos=item.z_pos,
                                                  x_rate=x_special_pos, y_rate=y_special_pos, z_rate=z_special_pos,
                                                  x_sigma=item.x_sigma, y_sigma=item.y_sigma, z_sigma=item.z_sigma,
                                                  first_epoch=item.first_epoch, last_epoch=item.last_epoch);
        special_pos_objects.append(newobj);
    return special_pos_objects;


def get_Helmert_parameters(xyz_velfieldA, xyz_velfieldB):
    """
    Find the 7-parameter Helmert transformation minimizing the misfit between two reference frames.
    The parameters convert SystemA into SystemB.
    xyz_velfieldA and xyz_velfieldB must be lists with identical stations.
    raw xyz position is in x_pos/y_pos/z_pos
    velocity-adjusted position is in x_rate/y_rate/z_rate.
    """

    def fun(v):
        vecA, vecB = [], [];
        [Cx, Cy, Cz, s, Rx, Ry, Rz] = [*v];
        s_mult = 1 + s*1e-6;
        for i in range(len(xyz_velfieldA)):
            newobs = np.array([xyz_velfieldB[i].x_rate, xyz_velfieldB[i].y_rate, xyz_velfieldB[i].z_rate]);
            transformed_Aobs_x = Cx + s_mult * (xyz_velfieldA[i].x_rate - Rz*xyz_velfieldA[i].y_rate +
                                                Ry*xyz_velfieldA[i].z_rate);
            transformed_Aobs_y = Cy + s_mult * (Rz*xyz_velfieldA[i].x_rate + xyz_velfieldA[i].y_rate -
                                                Rx*xyz_velfieldA[i].z_rate);
            transformed_Aobs_z = Cz + s_mult * (-Ry*xyz_velfieldA[i].x_rate + Rx*xyz_velfieldA[i].y_rate +
                                                xyz_velfieldA[i].z_rate);
            transformed_obsA = np.array([transformed_Aobs_x, transformed_Aobs_y, transformed_Aobs_z]);
            vecB = np.concatenate((vecB, newobs), axis=0);
            vecA = np.concatenate((vecA, transformed_obsA), axis=0);
        residuals = np.subtract(vecB, vecA);
        return residuals;

    Htrans_starting_guess = np.array([0, 0, 0, 1, 0, 0, 0]);
    # solve nonlinear least squares problem!
    response = scipy.optimize.least_squares(fun, Htrans_starting_guess, loss='soft_l1');
    Htrans_optimal = response.x;  # the 7 parameters of the Helmert transformation.
    print(Htrans_optimal);
    return Htrans_optimal;


def postproc_after_helmert(xyz_velfield):
    """
    Take a special xyz velfield used for Helmert transformations
    Return a normal velfield in ENU
    """
    multiplier = 100;
    enu_station_list = [];
    for item in xyz_velfield:
        res_vel_x = (item.x_rate - item.x_pos) / multiplier;  # converting back into position and velocity
        res_vel_y = (item.y_rate - item.y_pos) / multiplier;
        res_vel_z = (item.z_rate - item.z_pos) / multiplier;
        lonlat = geo_conv.xyz2llh(np.array([[item.x_pos, item.y_pos, item.z_pos], ]));
        llh_origin_simple = np.array([lonlat[0][0], lonlat[0][1], 0]);  # simpler numpy array
        xyz_vel = 1000 * np.array([[res_vel_x, res_vel_y, res_vel_z], ]);  # now in mm
        xyz_stds = 1000 * np.array([item.x_sigma, item.y_sigma, item.z_sigma]);  # now in mm
        ecov = np.diag(xyz_stds);  # a 3x3 np array matrix with covariances on the diagonal

        enu_vel, enu_cov = geo_conv.xyz2enu(xyz_vel, llh_origin_simple, ecov);  # covariances go here

        new_station_vel = gps_io_functions.Station_Vel(name=item.name, elon=lonlat[0][0], nlat=lonlat[0][1],
                                                       e=enu_vel[0][0], n=enu_vel[0][1], u=enu_vel[0][2],
                                                       se=enu_cov[0][0], sn=enu_cov[1][1], su=enu_cov[2][2],
                                                       first_epoch=0, last_epoch=0, refframe=0, proccenter=0,
                                                       subnetwork=0, survey=0, meas_type=None);
        enu_station_list.append(new_station_vel);
    return enu_station_list;


def get_bounding_box(velfield, border=0.1):
    """Get a bounding box for a list of station_vels.  Returns a list of [E, W, S, N]"""
    lons = [x.elon for x in velfield]
    lats = [x.nlat for x in velfield]
    bbox = [np.min(lons)-border, np.max(lons)+border, np.min(lats)-border, np.max(lats)+border];
    return bbox;


"""
Steps to solving the Helmert transformation:
DONE: Get a list of 551 SCEC_SYS_B/MIDAS_SYS_A common stations, with known lon/lat coordinates
DONE: Measure their initial ENU velocity misfits for good measure
DONE: Convert all coordinates to XYZ
DONE: Convert all velocities to XYZ
DOME: Add some multiple of each XYZ velocity (like 10 years, to avoid floating point math) to each XYZ position
DONE: Perform H_Trans estimator for SYS_A to SYS_B (scipy.optimize.least_squares, minimizing 551*3 residuals)
DONE: Return 7 Parameters
DONE: Apply H_Trans to MIDAS_SYS_A, converting it into MIDAS_SYS_B
DONE: Subtract the stations' XYZ coordinates from SYS_B results to get residual XYZ velocities in SYS_B
DONE: Convert residual velocities back into ENU
DONE: Divide by multiple of years used to get yearly velocities (like 10 years)
DONE: Measure the SCEC_SYS_B/MIDAS_SYS_B ENU velocity misfits, which are hopefully smaller than before H_Trans
ANS: RMS for all 551 stations went from 9.1 to 9.0 mm. 
"""
