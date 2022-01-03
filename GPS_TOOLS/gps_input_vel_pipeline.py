"""Functions to import velocity fields"""

import sys
from . import gps_io_functions


def import_velfield(gps_config_file, network='pbo', refframe='ITRF', sub_network=''):
    """Read a velocity field from a certain network and refframe"""
    myParams = gps_io_functions.read_config_file(gps_config_file);
    if network == 'pbo':
        pbo_velfile = get_pbo_velfile(myParams.pbo_velocities, refframe);
        [myVelocities] = gps_io_functions.read_pbo_vel_file(pbo_velfile);
    elif network == 'cwu':
        cwu_velfile = get_cwu_velfile(myParams.pbo_velocities, refframe);
        [myVelocities] = gps_io_functions.read_pbo_vel_file_format(cwu_velfile);
    elif network == 'unr':
        unr_velfile = get_unr_velfile(myParams.unr_velocities, refframe);
        [myVelocities] = gps_io_functions.read_unr_vel_file(unr_velfile, myParams.unr_coords_file);
    elif network[0:4] == 'usgs':
        if len(network) == 4:  # if the network is just usgs
            usgs_velfile = get_usgs_velfile(myParams.usgs_vel_dir, refframe, sub_network);
        else:  # if the network has format similar to 'usgs-Pacific_Northwest'
            sub_network = network.split('-')[1];
            usgs_velfile = get_usgs_velfile(myParams.usgs_vel_dir, refframe, sub_network);
        [myVelocities] = gps_io_functions.read_usgs_velfile(usgs_velfile, myParams.usgs_cache_file);
    else:
        print("Error! Invalid choice of network [pick one of pbo/cwu/unr/usgs]");
        sys.exit(0);
    return myVelocities;


def get_pbo_velfile(velocities_dir, refframe='NA'):
    if refframe == 'NA':
        velfile = velocities_dir + 'NAM08_pbovelfile_feb2018.txt';
    else:
        velfile = velocities_dir + 'IGS08_pbovelfile_feb2018.txt';
    return velfile;


def get_cwu_velfile(velocities_dir, refframe='NA'):
    if refframe == 'NA':
        velfile = velocities_dir + 'cwu.final_nam14.vel'
    else:
        velfile = velocities_dir + 'cwu.final_igs14.vel'
    return velfile;


def get_unr_velfile(velocities_dir, refframe):
    if refframe == 'NA':
        velfile = velocities_dir + 'midas.NA_nov2021.txt'
    else:
        velfile = velocities_dir + 'midas.IGS14_nov2021.txt'
    return velfile;


def get_usgs_velfile(velocities_dir, refframe, sub_network):
    if sub_network == '':
        print("Error! Must provide sub-network for USGS velocity field");
        sys.exit(0);
    if refframe == 'NA':
        velfile = velocities_dir + sub_network + '/NAM_' + sub_network + '_vels.txt';
    else:
        velfile = velocities_dir + sub_network + '/ITRF_' + sub_network + '_vels.txt';
    return velfile;
