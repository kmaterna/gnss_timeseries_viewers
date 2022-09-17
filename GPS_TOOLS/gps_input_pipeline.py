"""Reading inputs for various time series and velocity files"""

import subprocess, sys, os, glob
import datetime as dt
from . import gps_io_functions, offsets


# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su',
# 'EQtimes']);  # in mm

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
    station_names = remove_blacklist(data_config_file, station_names);
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
    Function to access the time-series reading library.
    refframe choices are NA and ITRF.
    """
    Params = gps_io_functions.read_config_file(data_config_file);
    filename, sub_network = pre_screen_datasource_paths(Params, station, datasource, refframe, sub_network);

    if datasource == 'unr':
        [myData, offset_obj, eq_obj] = get_unr(Params, filename);  # UNR data format
    elif datasource in ['pbo', 'cwu', 'nmt']:
        [myData, offset_obj, eq_obj] = get_pbo_type(Params, filename, station, datasource);  # PBO data format
    elif datasource == 'usgs':
        [myData, offset_obj, eq_obj] = get_usgs(Params, filename, station, refframe, sub_network);  # USGS data
    elif datasource in ['gldas', 'nldas', 'noah025']:
        [myData, offset_obj, eq_obj] = get_gldas(Params, filename);  # hydro
    elif datasource == 'grace':
        [myData, offset_obj, eq_obj] = get_grace(filename);  # GRACE model
    elif datasource == 'lsdm':
        [myData, offset_obj, eq_obj] = get_lsdm(filename);  # LSDM model
    else:
        print('Error! data source "%s" not recognized. Returning empty object. ' % datasource);
        return [[], [], []];  # Error code.
    return [myData, offset_obj, eq_obj];


def get_unr(Params, filename):
    [myData] = gps_io_functions.read_UNR_magnet_ts_file(filename, Params["unr"]["directory"]+Params["unr"]["coords_file"]);
    Offsets = get_unr_offsets(myData, Params["unr"]["directory"]+Params["unr"]["offsets_file"],
                              Params["unr"]["directory"]+Params["unr"]["user_offsets_file"]);
    Earthquakes = get_unr_earthquakes(myData, Params["unr"]["directory"]+Params["unr"]["offsets_file"],
                                      Params["unr"]["directory"]+Params["unr"]["user_offsets_file"]);
    return [myData, Offsets, Earthquakes];

def get_pbo_type(Params, pbo_filename, station, database):
    [myData] = gps_io_functions.read_pbo_pos_file(pbo_filename);  # PBO data format
    Offsets = get_pbo_offsets(station, Params[database]["directory"]+Params[database]["offsets_dir"]);
    Earthquakes = get_pbo_earthquakes(station, Params[database]["directory"]+Params[database]["earthquakes_dir"]);
    return [myData, Offsets, Earthquakes];

def get_usgs(Params, filename, station, refframe, sub_network):
    [myData] = gps_io_functions.read_USGS_ts_file(filename);
    usgs_offsets_dir = Params['usgs']['directory'] + Params['usgs']['offsets_dir'];
    Offsets = get_usgs_offsets(station, usgs_offsets_dir, sub_network, refframe, typekey='antenna/other');
    Earthquakes = get_usgs_offsets(station, usgs_offsets_dir, sub_network, refframe, typekey='earthquake');
    return [myData, Offsets, Earthquakes];

def get_gldas(Params, filename):
    coords_file = Params["unr"]["directory"]+Params["unr"]["coords_file"];
    [myData] = gps_io_functions.read_pbo_hydro_file(filename, coords_file);
    return [myData, [], []];

def get_grace(filename):
    [myData] = gps_io_functions.read_grace(filename);
    return [myData, [], []];

def get_lsdm(filename):
    [myData] = gps_io_functions.read_lsdm_file(filename);
    return [myData, [], []];


def pre_screen_datasource_paths(Params, station, inps='pbo', refframe="NA", sub_network=''):
    """ A defensive programming function that quits if it doesn't find the right file."""
    print("\nStation %s: " % station);

    if refframe not in ['NA', 'ITRF']:
        print("Error! Reference frame doesn't match available ones [NA, ITRF]. Choose again"); sys.exit(1);

    # CHECK IF STATION EXISTS IN JUST ONE SUBNETWORK FOR CONVENIENCE
    if inps == 'usgs' and sub_network == '':
        network_list = query_usgs_network_name(station, Params['usgs']['directory']+Params['usgs']['gps_ts_dir']);
        if len(network_list) == 1:
            sub_network = network_list[0].split('/')[-1];
        else:
            print("ERROR! User must select one sub-network for USGS time series. Exiting. ");
            sys.exit(1);

    # Path-setting through a look-up table
    if inps == 'unr':
        na_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.NA.tenv3';
        itrf_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.IGS14.tenv3';
    elif inps in ['cwu', 'pbo', 'nmt']:
        na_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.' + inps + '.final_nam14.pos';
        itrf_ts_filename = Params[inps]["directory"] + Params[inps]["gps_ts_dir"] + station + '.' + inps + '.final_igs14.pos';
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
            query_usgs_network_name(station, Params["usgs"]["directory"]);
        print("Exiting immediately...");
        sys.exit(1);
    return final_filename, sub_network;


def query_usgs_network_name(station, gps_ts_dir):
    """
    Given that USGS puts stations into networks, which network is a given station a member of?
    This function may print more than one.
    """
    network_list = [];
    directories = glob.glob(gps_ts_dir+"/*");
    print("Querying for station %s in USGS sub-networks... " % station);
    for item in directories:
        if os.path.isfile(item+'/'+station.lower()+'_NAfixed.rneu'):
            print("Found %s in %s" % (station, item) );
            network_list.append(item);
    print("")
    return network_list;


def remove_blacklist(data_config_file, stations):
    Params = gps_io_functions.read_config_file(data_config_file);
    new_stations = [];
    blacklisted_stations = gps_io_functions.read_blacklist(Params["blacklist"]);
    for station in stations:
        if not (station in blacklisted_stations):
            new_stations.append(station);
        else:
            print("Excluding station %s due to blacklist_file %s" % (station, Params["blacklist"]));
    return new_stations;


# THE GUTS --------------------------------------
def get_unr_offsets(Data0, offsets_file, user_offsets_file):
    """ grep -1 is the code for antenna and reference frame offsets"""
    print("Offset table for station %s:" % Data0.name);
    evdts1 = search_file_for_unr_offsets(Data0.name, offsets_file, mode=1);
    evdts2 = search_file_for_unr_offsets(Data0.name, user_offsets_file, mode=1);
    UNR_offsets = offsets.solve_for_offsets(Data0, evdts1 + evdts2);
    offsets.print_offset_object(UNR_offsets);
    return UNR_offsets;


def get_unr_earthquakes(Data0, offsets_file, user_offsets_file):
    """ grep -2 is the code for earthquake offsets """
    print("Earthquakes table for station %s:" % Data0.name);
    evdts1 = search_file_for_unr_offsets(Data0.name, offsets_file, mode=2);
    evdts2 = search_file_for_unr_offsets(Data0.name, user_offsets_file, mode=2);
    UNR_earthquakes = offsets.solve_for_offsets(Data0, evdts1 + evdts2);
    offsets.print_offset_object(UNR_earthquakes);
    return UNR_earthquakes;


def search_file_for_unr_offsets(station, offsets_file, mode=2):
    """
    Mode1 is antenna etc.
    Mode2 is earthquakes
    """
    try:
        table = subprocess.check_output(
            "grep -E '" + station + "  [0-9]{2}[A-Z]{3}[0-9]{2}  "+str(mode)+"' " + offsets_file, shell=True);
        table = table.decode();  # needed when switching to python 3
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    evdts = parse_table_unr(table);
    return evdts;


def get_pbo_offsets(station, offsets_dir):
    print("Offset table for station %s:" % station);
    PBO_offsets = [];
    table = search_pbo_offset_table(station, offsets_dir, "cwu*.off")  # Read the offset table
    [e_offsets, n_offsets, u_offsets, evdts] = parse_antenna_table_pbo(table);
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        PBO_offsets.append(offi);
    return PBO_offsets;


def get_pbo_earthquakes(station, earthquakes_dir):
    print("Earthquake table for station %s:" % station);
    PBO_earthquakes = [];
    table = search_pbo_offset_table(station, earthquakes_dir, "pbo*kalts.evt")  # Read the events table
    [e_offsets, n_offsets, u_offsets, evdts] = parse_earthquake_table_pbo(table);
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        PBO_earthquakes.append(offi);
    return PBO_earthquakes;


def get_cwu_earthquakes(station, earthquakes_dir):
    print("Earthquake table for station %s:" % station);
    CWU_earthquakes = [];
    table = search_pbo_offset_table(station, earthquakes_dir, "cwu*kalts.evt")  # Read the events table
    [e_offsets, n_offsets, u_offsets, evdts] = parse_earthquake_table_pbo(table);
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        CWU_earthquakes.append(offi);
    return CWU_earthquakes;


def search_pbo_offset_table(station, offset_dir, search_path):
    try:
        table = subprocess.check_output("grep " + station + " " + offset_dir + search_path, shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    print(table);
    return table;


def get_usgs_offsets(station, usgs_offsets_dir, sub_network, refframe, typekey='antenna/other'):
    # For earthquakes, antenna changes, and other offsets
    print(typekey, "table for station %s in %s :" % (station, refframe));
    USGS_offsets = [];
    offset_file = usgs_offsets_dir + sub_network + '/' + refframe + '_' + sub_network + '_offsets.txt';
    # Read the offset table
    try:
        table = subprocess.check_output("grep " + station + " " + offset_file, shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    [e_offsets, n_offsets, u_offsets, evdts] = parse_table_usgs(table, typekey);
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        USGS_offsets.append(offi);
    offsets.print_offset_object(USGS_offsets);
    return USGS_offsets;


# TABLE INPUTS ---------------------------
def parse_antenna_table_pbo(table):
    e_offsets, n_offsets, u_offsets, evdts = [], [], [], [];
    if len(table) == 0:
        return [e_offsets, n_offsets, u_offsets, evdts];
    table_rows = table.split('\n');
    for line in table_rows:
        if "EQ" in line:
            continue;
        else:
            print(line);
        if len(line) == 0:
            continue;  # if we're at the end, move on.
        words = line.split();
        yyyy, mm, dd = words[1], words[2], words[3];
        e_offsets.append(float(words[8]));  # in mm
        n_offsets.append(float(words[6]));
        u_offsets.append(float(words[10]));
        evdts.append(dt.datetime.strptime(yyyy + mm + dd, "%Y%m%d"));
    return [e_offsets, n_offsets, u_offsets, evdts];


def parse_earthquake_table_pbo(table):
    e_offsets, n_offsets, u_offsets, evdts = [], [], [], [];
    if len(table) == 0:
        return [e_offsets, n_offsets, u_offsets, evdts];
    tablesplit = table.split('\n');
    for item in tablesplit:  # for each earthquake
        if len(item) == 0:
            continue;  # if we're at the end, move on.
        words = item.split();
        filename = words[0];
        e_offsets.append(float(words[3]));  # in mm
        n_offsets.append(float(words[4]));
        u_offsets.append(float(words[8]));
        evdate = filename.split('/')[-1];
        evdate = evdate[4:10];
        year = evdate[0:2];
        month = evdate[2:4];
        day = evdate[4:6];
        year = "20" + year;
        evdts.append(dt.datetime.strptime(year + month + day, "%Y%m%d"));
    return [e_offsets, n_offsets, u_offsets, evdts];


def parse_table_usgs(table, offset_type):
    e_offsets, n_offsets, u_offsets, evdts = [], [], [], [];
    if len(table) == 0:
        return [e_offsets, n_offsets, u_offsets, evdts];
    tablesplit = table.split('\n');
    for item in tablesplit:  # for each earthquake
        if offset_type == 'earthquake':
            if 'earthquake' in item:
                words = item.split();
                if len(words) < 8:
                    continue;
                evdts.append(dt.datetime.strptime(words[1], "%Y-%m-%d"));
                n_offsets.append(float(words[3]))
                e_offsets.append(float(words[5]))
                u_offsets.append(float(words[7]))
        else:
            if 'earthquake' not in item:
                words = item.split();
                if len(words) < 8:
                    continue;
                evdts.append(dt.datetime.strptime(words[1], "%Y-%m-%d"));
                n_offsets.append(float(words[3]))
                e_offsets.append(float(words[5]))
                u_offsets.append(float(words[7]))

    return [e_offsets, n_offsets, u_offsets, evdts];


def parse_table_unr(table):
    """Here we extract all the antenna or earthquake offsets from the UNR table."""
    evdts = [];
    if len(table) == 0:  # if an empty list.
        return [];
    tablesplit = table.split('\n');
    for item in tablesplit:
        if len(item) == 0:
            continue;
        words = item.split();
        datestring = words[1];
        mydt = get_datetime_from_unrfile(datestring);
        if mydt not in evdts:  # we don't need redundant entries on the same date.
            # What if it happens within a week of each other?  Haven't figured this out yet.
            evdts.append(mydt);
    return evdts;


def get_datetime_from_unrfile(input_string):
    """Turns something like "12FEB13" into datetime.dt object for 2012-02-13"""
    year = input_string[0:2];
    if int(year) >= 80:
        year = '19' + year;
    else:
        year = '20' + year;
    mydt = dt.datetime.strptime(year + input_string[2:], "%Y%b%d");
    return mydt;
