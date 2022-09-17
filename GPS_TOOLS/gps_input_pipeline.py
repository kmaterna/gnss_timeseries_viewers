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
    datasource, sub_network = pre_screen_datasource(data_config_file, station, datasource, refframe, sub_network);

    if datasource == 'pbo':
        [myData, offset_obj, eq_obj] = get_pbo(data_config_file, station, refframe);  # PBO data format
    elif datasource == 'unr':
        [myData, offset_obj, eq_obj] = get_unr(data_config_file, station, refframe);  # UNR data format
    elif datasource == 'cwu':
        [myData, offset_obj, eq_obj] = get_cwu(data_config_file, station, refframe);  # CWU data
    elif datasource == 'nmt':
        [myData, offset_obj, eq_obj] = get_nmt(data_config_file, station, refframe);  # NMT data
    elif datasource == 'usgs':
        [myData, offset_obj, eq_obj] = get_usgs(data_config_file, station, sub_network, refframe);  # USGS data
    elif datasource == 'gldas':
        [myData, offset_obj, eq_obj] = get_gldas(data_config_file, station);  # GLDAS hydro
    elif datasource == 'nldas':
        [myData, offset_obj, eq_obj] = get_nldas(data_config_file, station);  # NLDAS hydro
    elif datasource == 'noah025':
        [myData, offset_obj, eq_obj] = get_noah025(data_config_file, station);  # NOAH0.25
    elif datasource == 'grace':
        [myData, offset_obj, eq_obj] = get_grace(data_config_file, station);  # GRACE model
    elif datasource == 'lsdm':
        [myData, offset_obj, eq_obj] = get_lsdm(data_config_file, station);  # LSDM model
    else:
        print('Error! data source "%s" not recognized. Returning empty object. ' % datasource);
        return [[], [], []];  # Error code.
    return [myData, offset_obj, eq_obj];


def get_unr(data_config_file, station, refframe="NA"):
    if refframe == "ITRF":
        reflabel = "IGS14";
    else:
        reflabel = "NA";
    Params = gps_io_functions.read_config_file(data_config_file);
    unr_filename = Params.unr_gps_dir + station + "." + reflabel + ".tenv3"
    [myData] = gps_io_functions.read_UNR_magnet_ts_file(unr_filename, Params.unr_coords_file);  # UNR data format
    Offsets = get_unr_offsets(myData, Params.unr_offsets_file, Params.unr_user_offsets_file);
    Earthquakes = get_unr_earthquakes(myData, Params.unr_offsets_file, Params.unr_user_offsets_file);
    return [myData, Offsets, Earthquakes];


def get_pbo(data_config_file, station, refframe="NA"):
    if refframe == "ITRF":
        reflabel = "igs08";
    else:
        reflabel = "nam08"
    Params = gps_io_functions.read_config_file(data_config_file);
    pbo_filename = Params.pbo_gps_dir + station + ".pbo.final_" + reflabel + ".pos"
    [myData] = gps_io_functions.read_pbo_pos_file(pbo_filename);  # PBO data format
    Offsets = get_pbo_offsets(station, Params.pbo_offsets_dir);
    Earthquakes = get_pbo_earthquakes(station, Params.pbo_earthquakes_dir);
    return [myData, Offsets, Earthquakes];


def get_cwu(data_config_file, station, refframe="NA"):
    if refframe == "ITRF":
        reflabel = "igs14";
    else:
        reflabel = "nam14"
    Params = gps_io_functions.read_config_file(data_config_file);
    pbo_filename = Params.pbo_gps_dir + station + ".cwu.final_" + reflabel + ".pos"
    [myData] = gps_io_functions.read_pbo_pos_file(pbo_filename);  # PBO data format
    Offsets = get_pbo_offsets(station, Params.pbo_offsets_dir);
    Earthquakes = get_cwu_earthquakes(station, Params.pbo_earthquakes_dir);
    return [myData, Offsets, Earthquakes];


def get_nmt(data_config_file, station, refframe="NA"):
    if refframe == "ITRF":
        reflabel = "igs08";
    else:
        reflabel = "nam08";
    Params = gps_io_functions.read_config_file(data_config_file);
    pbo_filename = Params.pbo_gps_dir + station + ".nmt.final_" + reflabel + ".pos"
    [myData] = gps_io_functions.read_pbo_pos_file(pbo_filename);  # PBO data format
    Offsets = get_pbo_offsets(station, Params.pbo_offsets_dir);
    Earthquakes = get_pbo_earthquakes(station, Params.pbo_earthquakes_dir);
    return [myData, Offsets, Earthquakes];


def get_usgs(data_config_file, station, sub_network, refframe="NA"):
    Params = gps_io_functions.read_config_file(data_config_file);
    if refframe == 'ITRF':
        reflabel = 'ITRF2008';
        reflabel_offsets = 'ITRF';
    else:
        reflabel = 'NAfixed';
        reflabel_offsets = 'NAM';
    usgs_filename = Params.usgs_gps_dir + sub_network + '/' + station.lower() + "_" + reflabel + ".rneu"
    [myData] = gps_io_functions.read_USGS_ts_file(usgs_filename);
    Offsets = get_usgs_offsets(station, Params.usgs_offsets_dir, sub_network, reflabel_offsets);
    Earthquakes = get_usgs_earthquakes(station, Params.usgs_offsets_dir, sub_network, reflabel_offsets);
    return [myData, Offsets, Earthquakes];


def get_gldas(data_config_file, station):
    Params = gps_io_functions.read_config_file(data_config_file);
    coords_file = Params.unr_coords_file;
    station_name_lower = station.lower();
    filename = Params.gldas_dir + station_name_lower + "_noah10_gldas2.hyd";
    [myData] = gps_io_functions.read_pbo_hydro_file(filename, coords_file);
    Offset = [];
    return [myData, Offset, Offset];


def get_nldas(data_config_file, station):
    Params = gps_io_functions.read_config_file(data_config_file);
    coords_file = Params.unr_coords_file;
    station_name_lower = station.lower();
    filename = Params.nldas_dir + station_name_lower + "_noah125_nldas2.hyd";
    [myData] = gps_io_functions.read_pbo_hydro_file(filename, coords_file);
    Offset = [];
    return [myData, Offset, Offset];


def get_noah025(data_config_file, station):
    Params = gps_io_functions.read_config_file(data_config_file);
    coords_file = Params.unr_coords_file;
    filename = Params.noah_dir + station + "_NOAH025.hyd";
    [myData] = gps_io_functions.read_pbo_hydro_file(filename, coords_file);
    Offset = [];
    return [myData, Offset, Offset];


def get_grace(data_config_file, station):
    Params = gps_io_functions.read_config_file(data_config_file);
    filename = Params.grace_dir + "scaled_" + station + "_PREM_model_ts.txt";
    [myData] = gps_io_functions.read_grace(filename);
    Offset = [];
    return [myData, Offset, Offset];


def get_lsdm(data_config_file, station):
    Params = gps_io_functions.read_config_file(data_config_file);
    filename = Params.lsdm_dir + station + "_LSDM_hydro.txt.txt";
    [myData] = gps_io_functions.read_lsdm_file(filename);
    Offset = [];
    return [myData, Offset, Offset];


def pre_screen_datasource(data_config_file, station, input_datasource='pbo', refframe="NA", sub_network=''):
    """ A defensive programming function that quits if it doesn't find the right file."""
    print("\nStation %s: " % station);
    Params = gps_io_functions.read_config_file(data_config_file);
    if refframe == "NA":
        unr_reflabel = "NA";
        pbo_reflabel = "nam08";
        cwu_reflabel = "nam14";
        usgs_reflabel = "NAfixed";
    elif refframe == "ITRF":
        unr_reflabel = "IGS14";
        pbo_reflabel = "igs08";
        cwu_reflabel = "igs14";
        usgs_reflabel = "ITRF2008";
    else:
        print("ERROR! Unrecognized reference frame (choices NA and ITRF)");
        sys.exit(1);

    if input_datasource == 'usgs' and sub_network == '':
        network_list = query_usgs_network_name(station, Params.usgs_gps_dir);
        if len(network_list) == 1:
            sub_network = network_list[0].split('/')[-1];
        else:
            print("ERROR! User must select one sub-network for USGS time series. Exiting. ");
            sys.exit(1);

    # Path setting
    if input_datasource == 'unr':
        filename = Params.unr_gps_dir + station + "." + unr_reflabel + ".tenv3";
    elif input_datasource == 'pbo':
        filename = Params.pbo_gps_dir + station + ".pbo.final_" + pbo_reflabel + ".pos";
    elif input_datasource == 'cwu':
        filename = Params.pbo_gps_dir + station + ".cwu.final_" + cwu_reflabel + ".pos";
    elif input_datasource == 'nmt':
        filename = Params.pbo_gps_dir + station + ".nmt.final_" + pbo_reflabel + ".pos";
    elif input_datasource == 'usgs':
        filename = Params.usgs_gps_dir + '/' + sub_network + '/' + station.lower() + "_" + usgs_reflabel + ".rneu";
    elif input_datasource == 'gldas':
        filename = Params.gldas_dir + station.lower() + "_noah10_gldas2.hyd";
    elif input_datasource == 'nldas':
        filename = Params.nldas_dir + station.lower() + "_noah125_nldas2.hyd";
    elif input_datasource == 'noah025':
        filename = Params.noah_dir + station + "_NOAH025.hyd";
    elif input_datasource == 'grace':
        filename = Params.grace_dir + "scaled_" + station + "_PREM_model_ts.txt";
    elif input_datasource == 'lsdm':
        filename = Params.lsdm_dir + station + "_LSDM_hydro.txt.txt";
    else:
        print("Error! Invalid input datasource %s" % input_datasource);
        sys.exit(1);

    # Determine if the file is found on the computer or not. Provide helpful suggestions if not.
    if os.path.isfile(filename):
        print("Found file %s from datasource %s %s " % (filename, input_datasource, sub_network));
    # If the file is not found on the computer:
    else:
        print("Error!  Cannot find %s in %s %s database." % (filename, input_datasource, sub_network) );
        if input_datasource == 'usgs':
            print("The station %s could be found in the following sub_networks instead: " % station);
            query_usgs_network_name(station, Params.usgs_gps_dir);
        print("Exiting immediately...");
        sys.exit(1);
    return input_datasource, sub_network;


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
    blacklisted_stations = gps_io_functions.read_blacklist(Params.blacklist);
    for station in stations:
        if not (station in blacklisted_stations):
            new_stations.append(station);
        else:
            print("Excluding station %s due to blacklist_file %s" % (station, Params.blacklist));
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
    try:
        table = subprocess.check_output("grep " + station + " " + offsets_dir + "cwu*.off", shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    print(table);
    [e_offsets, n_offsets, u_offsets, evdts] = parse_antenna_table_pbo(table);
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        PBO_offsets.append(offi);
    return PBO_offsets;


def get_pbo_earthquakes(station, earthquakes_dir):
    print("Earthquake table for station %s:" % station);
    PBO_earthquakes = [];
    # Read the offset table
    try:
        table = subprocess.check_output("grep " + station + " " + earthquakes_dir + "pbo*kalts.evt", shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    print(table);
    [e_offsets, n_offsets, u_offsets, evdts] = parse_earthquake_table_pbo(table);
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        PBO_earthquakes.append(offi);
    return PBO_earthquakes;


def get_cwu_earthquakes(station, earthquakes_dir):
    print("Earthquake table for station %s:" % station);
    CWU_earthquakes = [];
    # Read the offset table
    try:
        table = subprocess.check_output("grep " + station + " " + earthquakes_dir + "cwu*kalts.evt", shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    print(table);
    [e_offsets, n_offsets, u_offsets, evdts] = parse_earthquake_table_pbo(table);
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        CWU_earthquakes.append(offi);
    return CWU_earthquakes;


def get_usgs_offsets(station, usgs_offsets_dir, sub_network, refframe):
    # For antenna changes and other offsets
    print("Offset table for station %s in %s :" % (station, refframe));
    USGS_offsets = [];
    offset_file = usgs_offsets_dir + sub_network + '/' + refframe + '_' + sub_network + '_offsets.txt';
    # Read the offset table
    try:
        table = subprocess.check_output("grep " + station + " " + offset_file, shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    [e_offsets, n_offsets, u_offsets, evdts] = parse_table_usgs(table, 'antenna/other');
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        USGS_offsets.append(offi);
    offsets.print_offset_object(USGS_offsets);
    return USGS_offsets;


def get_usgs_earthquakes(station, usgs_offsets_dir, sub_network, refframe):
    # For earthquake offsets
    print("Earthquake table for station %s in %s :" % (station, refframe));
    USGS_earthquakes = [];
    offset_file = usgs_offsets_dir + sub_network + '/' + refframe + '_' + sub_network + '_offsets.txt';
    # Read the offset table
    try:
        table = subprocess.check_output("grep " + station + " " + offset_file, shell=True);
    except subprocess.CalledProcessError:  # if we have no earthquakes in the event files...
        table = [];
    if len(table) > 0:
        table = table.decode();  # needed when switching to python 3
    [e_offsets, n_offsets, u_offsets, evdts] = parse_table_usgs(table, 'earthquake');
    for i in range(len(e_offsets)):
        offi = offsets.Offsets(e_offsets=e_offsets[i], n_offsets=n_offsets[i], u_offsets=u_offsets[i], evdts=evdts[i]);
        USGS_earthquakes.append(offi);
    offsets.print_offset_object(USGS_earthquakes);
    return USGS_earthquakes;


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
        yyyy = words[1];
        mm = words[2];
        dd = words[3];
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
