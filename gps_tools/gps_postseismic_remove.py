"""
A function to remove postseismic deformation via existing model time series
The first model is from Hines et al., JGR, 2016
"""
from gps_tools.file_io import config_io, io_nota
import gps_tools.gps_objects
import numpy as np
import os
from . import offsets, gps_ts_functions


# HELPER FUNCTIONS #
def get_station_hines(station_name, data_config_file):
    """
    Reads a model time series from Hines et al., JGR, 2016 into a gps_object.
    The time series are stored on disk within a directory in the /contrib part of the GNSS data directory.
    """
    system_params = config_io.read_config_file(data_config_file);
    model_dir = "/Contrib_Data/Remove_postseismic/Hines/Stations/"
    model_file = system_params["gps_data_dir"] + model_dir + station_name + "_psmodel.pos";
    # This is stored in general_gps_dir because it's on my system, but may not be on general systems. 
    if os.path.isfile(model_file):
        [Data0] = io_nota.read_pbo_pos_file(model_file);
        return Data0;
    else:
        print("ERROR: Cannot remove postseismic because file does not exist; file %s" % model_file);
        return None;


def remove_by_model(data_obj, model_obj, starttime1, endtime1, starttime2, endtime2):
    """
    Remove a postseismic transient model from a gps_object,
    using a model time series formatted the same way as a gps_object.
    starttime and endtime parameters have to do with fixing the edges of the time spanned by the model.
    """
    if not model_obj:  # if None
        return data_obj;
    if len(data_obj.dtarray) == 0:  # if empty
        return data_obj;

    if model_obj.dtarray[-1] < data_obj.dtarray[-1]:
        print("\nWARNING! Trying to use a short postseismic model to fix a long GNSS time series.");
        print("PROBLEMS MAY OCCUR- tread carefully!!\n\n");

    # These will be the same size.
    Data0, model = gps_ts_functions.pair_gps_model_keeping_gps(data_obj, model_obj);
    # pair_gps_model_keeping_gps leaves data outside of the model timespan
    # Data0, model = gps_ts_functions.pair_gps_model(Data0, model_data);  # removes data outside of the model timespan.

    # Subtract model from data.
    dtarray = Data0.dtarray;
    dE_gps = np.array(np.subtract(Data0.dE, model.dE));
    dN_gps = np.array(np.subtract(Data0.dN, model.dN));
    dU_gps = np.array(np.subtract(Data0.dU, model.dU));

    # In this method, we correct for offsets at the beginning and end of the modeled time series.
    interval1 = [starttime1, endtime1];
    e_offset1 = offsets.fit_single_offset(dtarray, dE_gps, interval1, 20);
    n_offset1 = offsets.fit_single_offset(dtarray, dN_gps, interval1, 20);
    v_offset1 = offsets.fit_single_offset(dtarray, dU_gps, interval1, 20);
    offsets1 = gps_tools.gps_objects.Offsets(e_offsets=e_offset1, n_offsets=n_offset1, u_offsets=v_offset1,
                                             evdts=starttime1);
    interval2 = [starttime2, endtime2];
    e_offset2 = offsets.fit_single_offset(dtarray, dE_gps, interval2, 20);
    n_offset2 = offsets.fit_single_offset(dtarray, dN_gps, interval2, 20);
    v_offset2 = offsets.fit_single_offset(dtarray, dU_gps, interval2, 20);
    offsets2 = gps_tools.gps_objects.Offsets(e_offsets=e_offset2, n_offsets=n_offset2, u_offsets=v_offset2,
                                             evdts=starttime2);
    offsets_obj = [offsets1, offsets2];

    corrected_data = gps_tools.gps_objects.Timeseries(name=Data0.name, coords=Data0.coords, dtarray=dtarray, dE=dE_gps,
                                                      dN=dN_gps, dU=dU_gps, Se=Data0.Se, Sn=Data0.Sn, Su=Data0.Su,
                                                      EQtimes=Data0.EQtimes);
    corrected_data = offsets.remove_offsets(corrected_data, offsets_obj);
    return corrected_data;
