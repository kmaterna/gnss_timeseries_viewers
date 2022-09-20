"""Reading inputs for time series files"""

import sys, os
import datetime as dt
from .file_io import io_magnet_unr, io_nota, io_other, io_usgs, config_io
from . import offsets, utilities


# Timeseries = namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su', 'EQtimes']);  # in mm

# ----------------------------------------
#  MULTI STATION DRIVERS               ---
# ----------------------------------------

def multi_station_inputs(station_names, blacklist_user, proc_center, refframe, data_config_file, distances=None,
                         must_include=(None, None)):
    """
    Returns a list of objects for time series data, offsets, and earthquakes
    Also kicks out stations when necessary based on blacklist or timing criteria (through must_include parameter)
    Keeps the distances object as metadata in case you pass it through.
    must_include is a 2-vector of datetime objects that present a window that must be included in the time series.
    default is not to impose the must_include criterion.
    """
    dataobj_list, offsetobj_list, eqobj_list, distances_surviving = [], [], [], [];
    Params = config_io.read_config_file(data_config_file);
    station_names = utilities.remove_blacklist(station_names, io_other.read_blacklist(Params["blacklist"]));  # filter
    for i in range(len(station_names)):
        if station_names[i] in blacklist_user:
            continue;  # a second type of optional blacklist
        else:
            [myData, offset_obj, eq_obj] = get_station_data(station_names[i], proc_center, data_config_file, refframe);

            if must_include[0] is not None:
                if myData.dtarray[-1] < must_include[0] or myData.dtarray[0] > must_include[1]:
                    # kicking out the stations that end early or start late.
                    print("Excluding station %s due to not including data inside %s:%s" % (
                        myData.name, dt.datetime.strftime(must_include[0], "%Y-%m-%d"),
                        dt.datetime.strftime(must_include[1], "%Y-%m-%d")));
                    continue;
            if myData:
                dataobj_list.append(myData);
                offsetobj_list.append(offset_obj);
                eqobj_list.append(eq_obj);
                if distances is not None:  # if you want to sort distances too
                    distances_surviving.append(distances[i]);
    return [dataobj_list, offsetobj_list, eqobj_list, distances_surviving];


# ----------------------------------------
# DRIVERS, CONFIGURE, AND FILE MASHING ---
# ----------------------------------------

def get_station_data(station, datasource, data_config_file, refframe="NA", sub_network=''):
    """
    Function to access the time-series reading library. refframe choices are NA and ITRF.

    :param station: string, 4-characters
    :param datasource: string, like 'pbo' or 'unr'
    :param data_config_file: string, filepath
    :param refframe: string, default 'NA', choices usually 'NA' or 'ITRF'
    :param sub_network: string, sometimes required for USGS data
    :returns: [TimeSeries, list of Offsets for EQs, list of Offsets for antenna changes]
    """
    Params = config_io.read_config_file(data_config_file);
    filename, sub_network = pre_screen_datasource_paths(Params, station, datasource, refframe, sub_network);

    if datasource == 'unr':
        [myData, offset_obj, eq_obj] = get_unr(Params, filename);  # UNR data format
    elif datasource in ['pbo', 'nmt']:
        [myData, offset_obj, eq_obj] = get_pbo_type(Params, filename, station, datasource);  # PBO data format
    elif datasource == 'cwu':
        [myData, offset_obj, eq_obj] = get_cwu(Params, filename, station, datasource);  # CWU data format
    elif datasource == 'usgs':
        [myData, offset_obj, eq_obj] = get_usgs(Params, filename, station, refframe, sub_network);  # USGS data
    elif datasource in ['gldas', 'nldas', 'noah025']:
        [myData, offset_obj, eq_obj] = get_gldas(filename);  # hydro
    elif datasource == 'grace':
        [myData, offset_obj, eq_obj] = get_grace(filename);  # GRACE model
    elif datasource == 'lsdm':
        [myData, offset_obj, eq_obj] = get_lsdm(filename);  # LSDM model
    else:
        print('Error! data source "%s" not recognized. Returning empty object. ' % datasource);
        return [[], [], []];  # Error code.
    return [myData, offset_obj, eq_obj];


def get_unr(Params, filename):
    myData = io_magnet_unr.read_UNR_magnet_ts_file(filename, Params["unr"]["directory"] + Params["unr"]["coords_file"]);
    Offsets = get_unr_offsets(myData, Params["unr"]["directory"]+Params["unr"]["offsets_file"],
                              Params["unr"]["directory"]+Params["unr"]["user_offsets_file"]);
    Earthquakes = get_unr_earthquakes(myData, Params["unr"]["directory"]+Params["unr"]["offsets_file"],
                                      Params["unr"]["directory"]+Params["unr"]["user_offsets_file"]);
    return [myData, Offsets, Earthquakes];

def get_pbo_type(Params, pbo_filename, station, database):
    myData = io_nota.read_pbo_pos_file(pbo_filename);  # PBO data format
    Offsets = get_nota_offsets(station, Params[database]["directory"]+Params[database]["offsets_dir"],
                               file_pattern='pbo*.off');
    Earthquakes = get_nota_earthquakes(station, Params[database]["directory"]+Params[database]["earthquakes_dir"],
                                       file_pattern="pbo*kalts.evt");
    return [myData, Offsets, Earthquakes];

def get_cwu(Params, pbo_filename, station, database):
    myData = io_nota.read_pbo_pos_file(pbo_filename);  # PBO data format
    Offsets = get_nota_offsets(station, Params[database]["directory"]+Params[database]["offsets_dir"],
                               file_pattern='cwu*.off');
    Earthquakes = get_nota_earthquakes(station, Params[database]["directory"]+Params[database]["earthquakes_dir"],
                                       file_pattern="cwu*kalts.evt");
    return [myData, Offsets, Earthquakes];

def get_usgs(Params, filename, station, refframe, sub_network):
    myData = io_usgs.read_USGS_ts_file(filename);
    usgs_offsets_dir = Params['usgs']['directory'] + Params['usgs']['offsets_dir'];
    Offsets = get_usgs_offsets(station, usgs_offsets_dir, sub_network, refframe, typekey='antenna/other');
    Earthquakes = get_usgs_offsets(station, usgs_offsets_dir, sub_network, refframe, typekey='earthquake');
    return [myData, Offsets, Earthquakes];

def get_gldas(filename):
    myData = io_nota.read_pbo_hydro_file(filename);
    return [myData, [], []];

def get_grace(filename):
    myData = io_other.read_grace(filename);
    return [myData, [], []];

def get_lsdm(filename):
    myData = io_other.read_lsdm_file(filename);
    return [myData, [], []];


def pre_screen_datasource_paths(Params, station, inps='pbo', refframe="NA", sub_network=''):
    """ A defensive programming function that quits if it doesn't find the right file."""
    print("\nStation %s: " % station);

    if refframe not in ['NA', 'ITRF']:
        print("Error! Reference frame doesn't match available ones [NA, ITRF]. Choose again"); sys.exit(1);

    # CHECK IF STATION EXISTS IN JUST ONE SUBNETWORK FOR CONVENIENCE
    if inps == 'usgs' and sub_network == '':
        network_list = io_usgs.query_usgs_network_name(station, Params['usgs']['directory'] +
                                                       Params['usgs']['gps_ts_dir']);
        if len(network_list) == 1:
            sub_network = network_list[0].split('/')[-1];
        else:
            print("ERROR! User must select one sub-network for USGS time series. Exiting. ");
            sys.exit(1);

    # Path-setting through a look-up table
    if inps == 'unr':
        na_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.NA.tenv3';
        itrf_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.IGS14.tenv3';
    elif inps in ['cwu']:
        na_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.' + inps + '.final_nam14.pos';
        itrf_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.' + inps + '.final_igs14.pos';
    elif inps in ['pbo', 'nmt']:
        na_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.' + inps + '.final_nam08.pos';
        itrf_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.' + inps + '.final_igs08.pos';
    elif inps == 'usgs':  # needs subnetwork
        na_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + sub_network + '/' + station.lower()+"_NAfixed.rneu";
        itrf_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + sub_network + '/' + station.lower() + "_ITRF2008.rneu";
    elif inps == 'gldas':
        na_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_noah10_gldas2.hyd'
        itrf_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_noah10_gldas2.hyd'
    elif inps == 'nldas':
        na_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_noah125_nldas2.hyd'
        itrf_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_noah125_nldas2.hyd'
    elif inps == 'noah025':
        na_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_NOAH025.hyd'
        itrf_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_NOAH025.hyd'
    elif inps == 'grace':
        na_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + 'scaled_????_PREM_model_ts.txt'
        itrf_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + 'scaled_????_PREM_model_ts.txt'
    elif inps == 'lsdm':
        na_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_LSDM_hydro.txt.txt'
        itrf_ts_filename = Params["hydro"]["directory"] + Params["hydro"][inps+"_dir"] + station.lower() + '_LSDM_hydro.txt.txt'
    else:
        print("Error! Invalid input datasource %s" % inps);
        sys.exit(1);

    # Return the right filename based on selected reference frame
    if refframe == "NA":
        final_filename = na_ts_filename;
    elif refframe == "ITRF":
        final_filename = itrf_ts_filename;
    else:
        final_filename = itrf_ts_filename;

    # Determine if the file is found on the computer or not. Provide helpful suggestions if not.
    if os.path.isfile(final_filename):
        print("Found file %s from datasource %s %s " % (final_filename, inps, sub_network));
    # If the file is not found on the computer:
    else:
        print("Error!  Cannot find %s in %s %s database." % (final_filename, inps, sub_network) );
        if inps == 'usgs':
            print("The station %s could be found in the following sub_networks instead: " % station);
            io_usgs.query_usgs_network_name(station, Params["usgs"]["directory"]);
        print("Exiting immediately...");
        sys.exit(1);
    return final_filename, sub_network;


# THE GUTS --------------------------------------
def get_unr_offsets(Data0, offsets_file, user_offsets_file):
    """ grep -1 is the code for antenna and reference frame offsets"""
    print("Offset table for station %s:" % Data0.name);
    evdts1 = io_magnet_unr.search_file_for_unr_offsets(Data0.name, offsets_file, mode=1);
    evdts2 = io_magnet_unr.search_file_for_unr_offsets(Data0.name, user_offsets_file, mode=1);
    UNR_offsets = offsets.solve_for_offsets(Data0, evdts1 + evdts2);
    offsets.print_offset_object(UNR_offsets);
    return UNR_offsets;


def get_unr_earthquakes(Data0, offsets_file, user_offsets_file):
    """ grep -2 is the code for earthquake offsets """
    print("Earthquakes table for station %s:" % Data0.name);
    evdts1 = io_magnet_unr.search_file_for_unr_offsets(Data0.name, offsets_file, mode=2);
    evdts2 = io_magnet_unr.search_file_for_unr_offsets(Data0.name, user_offsets_file, mode=2);
    UNR_earthquakes = offsets.solve_for_offsets(Data0, evdts1 + evdts2);
    offsets.print_offset_object(UNR_earthquakes);
    return UNR_earthquakes;


def get_nota_offsets(station_name, offsets_dir, file_pattern):
    print("Offset table for station %s:" % station_name);
    table = io_nota.search_files_for_nota_offsets(station_name, offsets_dir, file_pattern)  # Read the offset table
    PBO_offsets = io_nota.parse_antenna_table_pbo(table);
    offsets.print_offset_object(PBO_offsets);
    return PBO_offsets;


def get_nota_earthquakes(station_name, earthquakes_dir, file_pattern):
    print("Earthquake table for station %s:" % station_name);
    table = io_nota.search_files_for_nota_offsets(station_name, earthquakes_dir, file_pattern)  # Read the events table
    PBO_earthquakes = io_nota.parse_earthquake_table_pbo(table);
    offsets.print_offset_object(PBO_earthquakes);
    return PBO_earthquakes;


def get_usgs_offsets(station, usgs_offsets_dir, sub_network, refframe, typekey='antenna/other'):
    # For earthquakes, antenna changes, and other offsets
    print(typekey, "table for station %s in %s :" % (station, refframe));
    offset_file = usgs_offsets_dir + sub_network + '/' + refframe + '_' + sub_network + '_offsets.txt';
    table = io_usgs.search_file_for_usgs_offsets(station, offset_file);
    USGS_offsets = io_usgs.parse_offset_table_usgs(table, typekey);
    offsets.print_offset_object(USGS_offsets);
    return USGS_offsets;
