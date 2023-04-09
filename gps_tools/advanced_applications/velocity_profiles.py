""" Profiles using pygmt """
from Tectonic_Utils.geodesy import haversine, insar_vector_functions
import numpy as np
import pygmt
import matplotlib.pyplot as plt
from .. import gps_objects as gps_objects
from .. import pygmt_plots


def project_onto_perpendicular_profile(velfield, startcoord, endcoord, width, vel_azimuth=None,
                                       fault_traces=None):
    """
    :param velfield: list of Station_Vel objects
    :param startcoord: tuple of beginning of profile (lon, lat)
    :param endcoord: tuple of end of profile (lon, lat)
    :param width: width of profile, in km
    :param vel_azimuth: optionally give azimuth (degrees east of north) of velocity component we want to profile
    :param fault_traces: List of fault trace files
    default is 90 degrees from azimuth of profile line itself.
    :returns: x-axis-distance, v_parallel, v_perp
    """

    profile_azimuth = haversine.calculate_initial_compass_bearing((startcoord[1], startcoord[0]),
                                                                  (endcoord[1], endcoord[0]));
    center_lon = np.nanmean([startcoord[0], endcoord[0]]);
    center_lat = np.nanmean([startcoord[1], endcoord[1]]);
    default_vel_azimuth = profile_azimuth - 90;
    print("profile azimuth: ", profile_azimuth, "degrees");
    if vel_azimuth is None:
        vel_azimuth = default_vel_azimuth;

    # make the pygmt projection
    # Build 783xN numpy array of stations and velocities
    lons = np.reshape(np.array([x.elon for x in velfield]), (len(velfield), 1));
    lats = np.reshape(np.array([x.nlat for x in velfield]), (len(velfield), 1));
    Ve = np.reshape(np.array([x.e for x in velfield]), (len(velfield), 1));
    Vn = np.reshape(np.array([x.n for x in velfield]), (len(velfield), 1));
    data = np.concatenate((lons, lats, Ve, Vn), axis=1);

    retval = pygmt.project(data, center=[center_lon, center_lat], azimuth=profile_azimuth,
                           width=[-width, width], unit=True);
    print(retval);

    # Pack up the selected stations in the profile
    selected_stations = [];
    for i in range(len(retval[0])):
        newstation = gps_objects.Station_Vel(name='', elon=retval[0][i], nlat=retval[1][i], e=retval[2][i],
                                             n=retval[3][i], u=0, proccenter=velfield[0].proccenter,
                                             refframe=velfield[0].refframe, subnetwork=velfield[0].subnetwork, se=0,
                                             sn=0, su=0, first_epoch=None, last_epoch=None, meas_type='gnss',
                                             survey=False);
        selected_stations.append(newstation);

    select_east = np.array(retval[2]);   # first column = east velocity
    select_north = np.array(retval[3]);   # second column = north velocity
    p = np.array(retval[4]);   # along-profile distance

    # Rotate velocities into azimuth
    rotation_theta = insar_vector_functions.bearing_to_cartesian(vel_azimuth);  # into deg north of east
    vel_parallel, vel_perp = [], [];
    for i in range(len(select_east)):
        xprime, yprime = insar_vector_functions.rotate_vector_by_angle(select_east[i], select_north[i], rotation_theta);
        vel_parallel.append(xprime);
        vel_perp.append(yprime);

    # Get the crossing locations of fault traces
    fault_x = [];
    if fault_traces:
        for item in fault_traces:
            fault_trace = np.loadtxt(item);
            retval = pygmt.project(data=fault_trace, center=[center_lon, center_lat], azimuth=profile_azimuth,
                                   width=[-15, 15], unit=True);
            fault_x.append(np.nanmean(retval[2]));   # append the average x-value of the crossings

    # Plotting stage: Map and Cartesian Plot.
    plot_profile_in_distance(p, vel_parallel, vel_perp, fault_xlocations=fault_x);
    pygmt_plots.map_velocity_profile(velfield, selected_stations, "profile.png", startcoord=startcoord,
                                     endcoord=endcoord, fault_traces=fault_traces);
    return;


# FAULT OBJECT: (x-loc, slip rate, locking depth)
def arctan_model(xmin, xmax, baselevel, fault_object_list):
    # FAULT OBJECT: (x-loc, slip rate, locking depth)
    x_array = np.arange(xmin, xmax, 0.1);
    y_array = baselevel + np.zeros(np.shape(x_array));
    for fault in fault_object_list:
        x0, sdot, d = fault[0], fault[1], fault[2]
        fault_contribution = sdot/2 - ((sdot/np.pi)*np.arctan((x_array-x0)/d));
        y_array = y_array + fault_contribution;
    return x_array, y_array;


def plot_profile_in_distance(x_distance, vel_parallel, vel_perp, fault_xlocations=()):
    """Make 1-D plot of velocities in the profile"""
    plt.figure(figsize=(10, 7), dpi=300);
    fontsize = 20;
    plt.plot(x_distance, vel_parallel, '.', linewidth=0, color='black', marker='s', label='Fault-parallel');
    plt.plot(x_distance, vel_perp, '.', linewidth=0, color='blue', marker='o', label='Fault-perpendicular');
    plt.xlabel("Distance along profile (km)", fontsize=fontsize);
    plt.ylabel('Velocity (mm/yr)', fontsize=fontsize);
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    if len(fault_xlocations) > 0:
        for item in fault_xlocations:
            top, bottom = [-10, 65];
            plt.plot([item, item], [top, bottom], '--', color='black');  # fault edges
    plt.ylim([-10, 59.9]);
    plt.xlim([-150, 251]);

    # SAF experiment:
    # Hard-coded parameters for a particular experiment on the Northern SAF
    saf_model0 = (fault_xlocations[0], 24, 12);
    maa_model0 = (fault_xlocations[1], 12, 12);
    bsf_model0 = (fault_xlocations[2], 6, 12);
    low_slip_rate_model = [saf_model0, maa_model0, bsf_model0]

    saf_model1 = (fault_xlocations[0], 15, 20);
    maa_model1 = (fault_xlocations[1], 22, 20);
    bsf_model1 = (fault_xlocations[2], 12, 20);
    high_slip_rate_model = [saf_model1, maa_model1, bsf_model1]

    xs, ys = arctan_model(-150, 250, 8, low_slip_rate_model);   # low slip rate case
    plt.plot(xs, ys, color='turquoise', linewidth=2, label='Elastic slip rates (24, 12, 6)');

    xs, ys = arctan_model(-150, 250, 8, high_slip_rate_model);   # high slip rate case
    plt.plot(xs, ys, color='purple', linewidth=2, label='VE slip rates (15, 22, 12)');

    plt.legend(fontsize=16);
    plt.title("Elastic half-space (arctangent) modeling, Point Arena profile", fontsize=20);
    plt.savefig("profile_in_distance.png");
    return;
