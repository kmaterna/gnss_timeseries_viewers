# A new function to remove postseismic deformation
# February 17, 2020
# The first model is from Hines et al., JGR, 2016

import datetime as dt
import os
import gps_io_functions
import offsets
import gps_ts_functions


# HELPER FUNCTIONS #
def get_station_hines(station_name, data_config_file):
    system_params = gps_io_functions.read_config_file(data_config_file);
    model_dir = "/Contrib_Data/Remove_postseismic/Hines/Stations/"
    model_file = system_params.general_gps_dir + model_dir + station_name + "_psmodel.pos";
    # This is stored in general_gps_dir because it's on my system, but may not be on general systems. 
    if os.path.isfile(model_file):
        [Data0] = gps_io_functions.read_pbo_pos_file(model_file);
        return Data0;
    else:
        print("ERROR: Cannot remove postseismic because file does not exist; file %s" % model_file);
        return None;


def remove_by_model(Data0, data_config_file):
    # Right now configured for the Hines data.
    starttime1 = dt.datetime.strptime("20100403", "%Y%m%d");
    endtime1 = dt.datetime.strptime("20100405", "%Y%m%d");
    starttime2 = dt.datetime.strptime("20200328", "%Y%m%d");
    endtime2 = dt.datetime.strptime("20200330", "%Y%m%d");

    # Input Hines data.
    model_data = get_station_hines(Data0.name, data_config_file);

    if not model_data:  # if None
        return Data0;

    if model_data.dtarray[-1] < Data0.dtarray[-1]:
        print("\nWARNING! Trying to use a short postseismic model to fix a long GNSS time series.");
        print("PROBLEMS MAY OCCUR- tread carefully!!\n\n");

    # These will be the same size.
    Data0, model = gps_ts_functions.pair_gps_model_keeping_gps(Data0, model_data);
    # pair_gps_model_keeping_gps leaves data outside of the model timespan
    # Data0, model = gps_ts_functions.pair_gps_model(Data0, model_data);  # removes data outside of the model timespan.

    # Subtract the model from the data.
    dtarray, dE_gps, dN_gps, dU_gps = [], [], [], [];
    Se_gps, Sn_gps, Su_gps = [], [], [];
    for i in range(len(Data0.dtarray)):
        dtarray.append(Data0.dtarray[i]);
        dE_gps.append(Data0.dE[i] - model.dE[i]);
        dN_gps.append(Data0.dN[i] - model.dN[i]);
        dU_gps.append(Data0.dU[i] - model.dU[i]);
        Se_gps.append(Data0.Se[i]);
        Sn_gps.append(Data0.Sn[i]);
        Su_gps.append(Data0.Su[i]);

    # In this method, we correct for offsets at the beginning and end of the modeled time series.
    interval1 = [starttime1, endtime1];
    east_offset1 = offsets.fit_single_offset(dtarray, dE_gps, interval1, 20);
    north_offset1 = offsets.fit_single_offset(dtarray, dN_gps, interval1, 20);
    vert_offset1 = offsets.fit_single_offset(dtarray, dU_gps, interval1, 20);
    interval2 = [starttime2, endtime2];
    east_offset2 = offsets.fit_single_offset(dtarray, dE_gps, interval2, 20);
    north_offset2 = offsets.fit_single_offset(dtarray, dN_gps, interval2, 20);
    vert_offset2 = offsets.fit_single_offset(dtarray, dU_gps, interval2, 20);

    offsets_obj = offsets.Offsets(e_offsets=[east_offset1, east_offset2],
                                  n_offsets=[north_offset1, north_offset2], u_offsets=[vert_offset1, vert_offset2],
                                  evdts=[starttime1, starttime2]);

    corrected_data = gps_io_functions.Timeseries(name=Data0.name, coords=Data0.coords, dtarray=dtarray, dE=dE_gps,
                                                 dN=dN_gps, dU=dU_gps, Se=Se_gps, Sn=Sn_gps, Su=Su_gps,
                                                 EQtimes=Data0.EQtimes);
    corrected_data = offsets.remove_offsets(corrected_data, offsets_obj);
    return corrected_data;
