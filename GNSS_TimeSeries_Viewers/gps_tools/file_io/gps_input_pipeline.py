"""Reading inputs for time series and their associated offset objects"""

from . import io_magnet_unr, io_nota, io_usgs
from .. import offsets


# ----------------------------------------
# DRIVERS, CONFIGURE, AND FILE MASHING ---
# ----------------------------------------

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
    all_evdts = evdts1 + evdts2;
    UNR_earthquakes = offsets.solve_for_offsets(Data0, all_evdts);
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
