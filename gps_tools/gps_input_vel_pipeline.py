"""Functions to import velocity fields"""

import sys
from gps_tools.file_io import io_magnet_unr, io_nota, io_usgs, config_io


def import_velfield(gps_config_file, network='pbo', refframe='ITRF', sub_network=''):
    """Read a velocity field from a certain network and refframe"""
    Params = config_io.read_config_file(gps_config_file);
    lookup_dict = build_lookup_dictionary(Params, network, sub_network);   # build a tranche of filenames
    if network == 'pbo':
        myVelocities = io_nota.read_pbo_vel_file(lookup_dict["pbo_" + refframe]);
    elif network == 'cwu':
        myVelocities = io_nota.read_pbo_vel_file_format(lookup_dict["cwu_" + refframe]);
    elif network == 'unr':
        myVelocities = io_magnet_unr.read_unr_vel_file(lookup_dict["unr_" + refframe], Params["unr"]["directory"] +
                                                       Params["unr"]["coords_file"]);
    elif network == 'usgs':
        myVelocities = io_usgs.read_usgs_velfile(lookup_dict["usgs_" + refframe], Params["usgs"]["directory"] +
                                                 Params["usgs"]["cache_file"]);
    else:
        print("Error! Invalid choice of network [pick one of pbo/cwu/unr/usgs]");
        sys.exit(0);
    return myVelocities;


def build_lookup_dictionary(Params, network, sub_network=''):
    """ Construct a dictionary of filepaths to velocity files. This can potentially grow if we add more frames."""
    lookup_dict = {};
    if network == 'pbo':
        lookup_dict["pbo_NA"] = Params["pbo"]["directory"] + Params["pbo"]["velocities_nam"];
        lookup_dict["pbo_ITRF"] = Params["pbo"]["directory"] + Params["pbo"]["velocities_itrf"];
    if network == 'cwu':
        lookup_dict["cwu_NA"] = Params["cwu"]["directory"] + Params["cwu"]["velocities_nam"];
        lookup_dict["cwu_ITRF"] = Params["cwu"]["directory"] + Params["cwu"]["velocities_itrf"];
    if network == 'unr':
        lookup_dict["unr_NA"] = Params["unr"]["directory"] + Params["unr"]["velocities_nam"];
        lookup_dict["unr_ITRF"] = Params["unr"]["directory"] + Params["unr"]["velocities_itrf"];
    if network == 'usgs':
        lookup_dict["usgs_NA"] = get_usgs_velfile(Params, 'NA', sub_network);
        lookup_dict["usgs_ITRF"] = get_usgs_velfile(Params, 'ITRF', sub_network);
    return lookup_dict;


def get_usgs_velfile(MyParams, refframe, sub_network):
    if sub_network == '':
        print("Error! Must provide sub-network for USGS velocity field");
        sys.exit(0);
    if refframe == 'NA':
        velfile = MyParams["usgs"]['directory'] + MyParams["usgs"]["vel_dir"] + sub_network + '/NAM_' + \
                  sub_network + '_vels.txt';
    else:
        velfile = MyParams["usgs"]['directory'] + MyParams["usgs"]["vel_dir"] + sub_network + '/ITRF_' + \
                  sub_network + '_vels.txt';
    return velfile;
