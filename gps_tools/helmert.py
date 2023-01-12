import numpy as np
import scipy.optimize
import Tectonic_Utils.geodesy.xyz2llh as geo_conv
from . import gps_objects, vel_functions


def create_Helmert_and_apply(sys_A_common_vels, sys_B_common_vels):
    """
    A full implementation and workflow of a 7-param Helmert transformation between two sets of common stations
    in two different reference frames.
    The result will be broadcast back in system A.

    :param sys_A_common_vels: list of station_vels in system A
    :param sys_B_common_vels: list of station_vels in system B
    :return: 7 Helmert params,
    """
    if len(sys_A_common_vels) != len(sys_B_common_vels):
        raise ValueError("Error! You are estimating a Helmert transformation between two non-matching vectors.");
    # Preparing for Estimating Helmert Transformation
    print("BEFORE HELMERT:", vel_functions.velocity_misfit_function(sys_A_common_vels, sys_B_common_vels),
          "mm RMS misfit");
    sysA_posform_xyz = prepare_velocities_for_helmert_trans(sys_A_common_vels);
    sysB_posform_xyz = prepare_velocities_for_helmert_trans(sys_B_common_vels);
    # ESTIMATING/APPLYING THE HELMERT TRANSFORMATION
    Hparams = get_Helmert_parameters(sysB_posform_xyz, sysA_posform_xyz);
    sysB_postH = Apply_Helmert_Transformation(sysB_posform_xyz, Hparams);
    # converting back to ENU
    B_enu_in_A = postproc_after_helmert(sysB_postH);
    A_enu_in_A = postproc_after_helmert(sysA_posform_xyz);
    print("AFTER HELMERT:", vel_functions.velocity_misfit_function(A_enu_in_A, B_enu_in_A), "mm RMS misfit");
    return Hparams, A_enu_in_A, B_enu_in_A;


def prepare_velocities_for_helmert_trans(velfield, multiplier=100):
    """
    Take regular Station_Vel objects (ENU).
    Then, we are using the "Station_Vel_XYZ" object in a clever way.
    x_pos, y_pos, z_pos are related to ONLY the lon/lat position
    x_rate, y_rate, z_rate are the position when a few years of velocity has been added
    x_sigma, y_sigma, z_sigma remain unchanged.
    This is going to be used for helmert transformation of positions into a new reference frame.
    Multiplier is the number of years for velocities, to avoid floating point and rounding errors
    """
    velfield_xyz = vel_functions.convert_enu_velfield_to_xyz(velfield);  # convert to xyz
    special_pos_objects = [];
    for item in velfield_xyz:
        x_special_pos = item.x_pos + multiplier * item.x_rate;
        y_special_pos = item.y_pos + multiplier * item.y_rate;
        z_special_pos = item.z_pos + multiplier * item.z_rate;
        newobj = gps_objects.Station_Vel_XYZ(name=item.name, x_pos=item.x_pos, y_pos=item.y_pos, z_pos=item.z_pos,
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
        s_mult = 1 + s * 1e-6;
        for i in range(len(xyz_velfieldA)):
            newobs = np.array([xyz_velfieldB[i].x_rate, xyz_velfieldB[i].y_rate, xyz_velfieldB[i].z_rate]);
            transformed_Aobs_x = Cx + s_mult * (xyz_velfieldA[i].x_rate - Rz * xyz_velfieldA[i].y_rate +
                                                Ry * xyz_velfieldA[i].z_rate);
            transformed_Aobs_y = Cy + s_mult * (Rz * xyz_velfieldA[i].x_rate + xyz_velfieldA[i].y_rate -
                                                Rx * xyz_velfieldA[i].z_rate);
            transformed_Aobs_z = Cz + s_mult * (-Ry * xyz_velfieldA[i].x_rate + Rx * xyz_velfieldA[i].y_rate +
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


def postproc_after_helmert(xyz_velfield, multiplier=100):
    """
    Take a special xyz velfield used for Helmert transformations
    Multiplier is the number of years for velocities, to avoid floating point and rounding errors
    Return a normal velfield in ENU
    """
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

        new_station_vel = gps_objects.Station_Vel(name=item.name, elon=lonlat[0][0], nlat=lonlat[0][1], e=enu_vel[0][0],
                                                  n=enu_vel[0][1], u=enu_vel[0][2], se=enu_cov[0][0], sn=enu_cov[1][1],
                                                  su=enu_cov[2][2], first_epoch=0, last_epoch=0, refframe=0,
                                                  proccenter=0, subnetwork=0, survey=0, meas_type=None);
        enu_station_list.append(new_station_vel);
    return enu_station_list;


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
    s_mult = 1 + s * 1e-6;
    H_Matrix = np.array([[1, -Rz, Ry],
                         [Rz, 1, -Rx],
                         [-Ry, Rx, 1]]);
    scaled_H_Matrix = np.multiply(H_Matrix, s_mult);
    translation_vector = np.array([[Cx, Cy, Cz]]);
    for item in xyz_velfieldA:
        # We store raw coordinates in x_pos, and velocity-considered coordinates in x_rate
        position_in_A = np.array([[item.x_rate, item.y_rate, item.z_rate]]);  # We use rate for transformed position
        position_in_B = np.dot(scaled_H_Matrix, position_in_A.T);  # Performing rotation/scaling transformation
        position_in_B = np.add(position_in_B, translation_vector.T);  # adding translation transformation
        position_in_B = position_in_B.T;
        xyz_stds = np.array([item.x_sigma, item.y_sigma, item.z_sigma]);
        ecov = np.diag(xyz_stds);  # a 3x3 np array matrix with covariances on the diagonal
        sigma_B = np.dot(np.dot(scaled_H_Matrix, ecov), scaled_H_Matrix.T);  # Transform sigmas
        sigma_B = np.multiply(sigma_B, np.square(s_mult));
        newobj = gps_objects.Station_Vel_XYZ(name=item.name, x_pos=item.x_pos, y_pos=item.y_pos, z_pos=item.z_pos,
                                             x_rate=position_in_B[0][0], y_rate=position_in_B[0][1],
                                             z_rate=position_in_B[0][2], x_sigma=sigma_B[0][0], y_sigma=sigma_B[1][1],
                                             z_sigma=sigma_B[2][2],
                                             first_epoch=item.first_epoch, last_epoch=item.last_epoch);
        xyz_velfieldB.append(newobj);
    return xyz_velfieldB;


def write_helmert_params(Hparams, filename):
    """
    :param Hparams: list of floats, seven params, [Cx, Cy, Cz, s, Rx, Ry, Rz]
    :param filename: string, filename
    """
    print("Writing %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("Cx: %e\n" % Hparams[0]);
    ofile.write("Cy: %e\n" % Hparams[1]);
    ofile.write("Cz: %e\n" % Hparams[2]);
    ofile.write("s: %e\n" % Hparams[3]);
    ofile.write("Rx: %e\n" % Hparams[4]);
    ofile.write("Ry: %e\n" % Hparams[5]);
    ofile.write("Rz: %e\n" % Hparams[6]);
    ofile.close();
    return;


"""
Steps to solving the Helmert transformation:
DONE: Get a list of SYS_B/SYS_A common stations, with known lon/lat coordinates
DONE: Measure their initial ENU velocity misfits for good measure
DONE: Convert all coordinates to XYZ
DONE: Convert all velocities to XYZ
DOME: Add some multiple of each XYZ velocity (like 100 years, to avoid floating point math) to each XYZ position
DONE: Perform H_Trans estimator for SYS_A to SYS_B (scipy.optimize.least_squares, minimizing N*3 residuals)
DONE: Return 7 Parameters
DONE: Apply H_Trans to MIDAS_SYS_A, converting it into MIDAS_SYS_B
DONE: Subtract the stations' XYZ coordinates from SYS_B results to get residual XYZ velocities in SYS_B
DONE: Convert residual velocities back into ENU
DONE: Divide by multiple of years used to get yearly velocities (like 100 years)
DONE: Measure the SYS_B/SYS_B ENU velocity misfits, which are hopefully smaller than before H_Trans
"""