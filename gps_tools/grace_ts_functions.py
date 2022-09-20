# GRACE FUNCTIONS
import numpy as np
import matplotlib.pyplot as plt
from . import utilities, gps_objects
from .file_io import io_other as io_other


def pair_GPSGRACE(GPS_TS, GRACE_TS):
    # This resamples the GRACE data to match GPS that is within the range of GRACE, and forms a common time axis.
    gps_decyear = utilities.get_float_times(GPS_TS.dtarray)
    grace_decyear = utilities.get_float_times(GRACE_TS.dtarray);  # the decimal years of all the grace obs points
    decyear, dtarray = [], [];
    north_gps, east_gps, vert_gps = [], [], [];
    N_err, E_err, V_err = [], [], [];
    for i in range(len(GPS_TS.dtarray)):  # this if-statement is happening because GPS is more current than GRACE
        if min(GRACE_TS.dtarray) < GPS_TS.dtarray[i] < max(GRACE_TS.dtarray):
            decyear.append(gps_decyear[i]);
            dtarray.append(GPS_TS.dtarray[i])
            north_gps.append(GPS_TS.dN[i]);
            east_gps.append(GPS_TS.dE[i]);
            vert_gps.append(GPS_TS.dU[i]);
            N_err.append(GPS_TS.Sn[i]);
            E_err.append(GPS_TS.Se[i]);
            V_err.append(GPS_TS.Su[i]);
    grace_u = np.interp(decyear, grace_decyear, GRACE_TS.dE);
    grace_v = np.interp(decyear, grace_decyear, GRACE_TS.dN);
    grace_w = np.interp(decyear, grace_decyear, GRACE_TS.dU);
    my_paired_ts = gps_objects.Paired_TS(dtarray=dtarray, north=north_gps, east=east_gps, vert=vert_gps, N_err=N_err,
                                         E_err=E_err, V_err=V_err, u=grace_u, v=grace_v, w=grace_w);
    return my_paired_ts;


def plot_grace(station_name, filename, out_dir):
    grace_ts = io_other.read_grace(filename);
    plt.figure();
    plt.plot_date(grace_ts.dtarray, grace_ts.u, '-b');
    plt.plot_date(grace_ts.dtarray, grace_ts.v, '-g');
    plt.plot_date(grace_ts.dtarray, grace_ts.w, '-r');
    plt.legend(['east', 'north', 'vertical']);
    plt.grid(True);
    plt.xlabel('Time');
    plt.ylabel('Displacement (mm)');
    plt.savefig(out_dir + station_name + "_gracets.eps");
    return;
