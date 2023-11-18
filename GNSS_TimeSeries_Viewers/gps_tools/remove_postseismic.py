"""
A function to remove postseismic deformation via existing model time series
"""
import numpy as np
from . import offsets, gps_ts_functions


def remove_by_model(data_obj, model_obj, starttime1, endtime1, starttime2, endtime2):
    """
    Remove a postseismic transient model from a gps_object,
    using a model time series formatted the same way as a gps_object.
    starttime and endtime parameters fix the edges of the time spanned by the model.
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
    offsets1 = offsets.Offset(e_offset=e_offset1, n_offset=n_offset1, u_offset=v_offset1, evdt=starttime1);
    interval2 = [starttime2, endtime2];
    e_offset2 = offsets.fit_single_offset(dtarray, dE_gps, interval2, 20);
    n_offset2 = offsets.fit_single_offset(dtarray, dN_gps, interval2, 20);
    v_offset2 = offsets.fit_single_offset(dtarray, dU_gps, interval2, 20);
    offsets2 = offsets.Offset(e_offset=e_offset2, n_offset=n_offset2, u_offset=v_offset2, evdt=starttime2);
    offsets_obj = [offsets1, offsets2];

    corrected_data = gps_ts_functions.Timeseries(name=Data0.name, coords=Data0.coords, dtarray=dtarray, dE=dE_gps,
                                                 dN=dN_gps, dU=dU_gps, Se=Data0.Se, Sn=Data0.Sn, Su=Data0.Su,
                                                 EQtimes=Data0.EQtimes);
    corrected_data = offsets.remove_offsets(corrected_data, offsets_obj);
    return corrected_data;
