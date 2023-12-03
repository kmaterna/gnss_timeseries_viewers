"""A set of functions that detrend or remove seasonals from TimeSeries objects, returning other TimeSeries objects

Non-Inclusive List of seasonal removal options:
   lssq: fits seasonals and linear trend by least squares inversion.
  notch: removes the 1-year and 6-month components by notch filter.
  grace: uses GRACE loading model interpolated between monthly points where available, and linear inversion where not.
    stl: uses a pre-computed look-up table for STL time series.
"""
import numpy as np
import datetime as dt
import glob, os, sys, subprocess
from . import gps_ts_functions, notch_filter, grace_ts_functions, utilities
from .gps_ts_functions import Timeseries
from .file_io import io_nota, io_other, config_io


def make_detrended_ts(Data, seasonals_remove, seasonals_type, data_config_file, remove_trend=1):
    """
    Generate a detrended and/or seasonally-removed time series. Seasonal fitting and de-trending in one function.
    There are several options for the seasonal removal (least-squares, notch filter, grace, etc.)
    """

    Params = config_io.read_config_file(data_config_file)

    if seasonals_remove == 0:
        print("Not removing seasonals.")
        [east_vel, north_vel, up_vel, _, _, _] = Data.get_slope(missing_fraction=0.6)
        east_params = [east_vel, 0, 0, 0, 0]   # Fit params definition: slope, a2(cos), a1(sin), s2, s1.
        north_params = [north_vel, 0, 0, 0, 0]
        up_params = [up_vel, 0, 0, 0, 0]
        trend_out = Data.detrend_data_by_value(east_params, north_params, up_params)
        trend_in = Data

    else:  # Going into different forms of seasonal removal.
        print("Removing seasonals by %s method." % seasonals_type)
        if seasonals_type == 'lssq':
            trend_out, trend_in = remove_seasonals_by_lssq(Data)

        elif seasonals_type == 'notch':
            trend_out, trend_in = remove_seasonals_by_notch(Data)

        elif seasonals_type == 'grace':
            trend_out, trend_in = remove_seasonals_by_GRACE(Data, Params['hydro']['grace_dir'])

        elif seasonals_type == 'stl':
            trend_out, trend_in = remove_seasonals_by_STL(Data, Params['hydro']['stl_dir'])

        elif seasonals_type == 'nldas':
            trend_out, trend_in = remove_seasonals_by_hydro(Data, Params['hydro']['nldas_dir'])

        elif seasonals_type == 'nldas_scaled':
            trend_out, trend_in = remove_seasonals_by_hydro(Data, Params['hydro']['nldas_dir'], scaling=True)

        elif seasonals_type == 'gldas':
            trend_out, trend_in = remove_seasonals_by_hydro(Data, Params['hydro']['gldas_dir'])

        elif seasonals_type == 'lsdm':
            trend_out, trend_in = remove_seasonals_by_german_load(Data, Params['hydro']['lsdm_dir'])

        elif seasonals_type == 'oroville':
            trend_out, trend_in = remove_seasonals_by_lakes(Data, Params['hydro']['lakes_dir'], 'oroville')

        elif seasonals_type == 'shasta':
            trend_out, trend_in = remove_seasonals_by_lakes(Data, Params['hydro']['lakes_dir'], 'shasta')

        else:
            print("Error: %s not supported as a seasonal removal type" % seasonals_type)
            print("The supported types are: lssq, grace, lsdm, nldas, nldas_scaled, gldas, notch, and stl")
            print("Exiting!\n")
            sys.exit(1)

    if remove_trend == 0:
        return trend_in
    else:
        return trend_out


def remove_seasonals_by_lssq(Data):
    [east_params, north_params, up_params] = Data.get_linear_annual_semiannual()
    trend_out = Data.detrend_data_by_value(east_params, north_params, up_params)
    trend_in = Data.remove_seasonal_by_value(east_params, north_params, up_params)
    return trend_out, trend_in


def simple_detrend_ts(dtarray, ts_array):
    """
    Remove a simple linear model from a time series

    :param dtarray: list of dts
    :param ts_array: 1d array of floats
    :return: 1d array of floats
    """
    decyear = utilities.get_float_times(dtarray)
    array_detrended = np.zeros(np.shape(decyear))
    model_coef = np.polyfit(decyear, ts_array, 1)[0]
    for i in range(len(ts_array)):
        array_detrended[i] = ts_array[i] - model_coef * decyear[i] - (ts_array[0] - model_coef * decyear[0])
    return array_detrended


def remove_seasonals_by_notch(Data):
    """
    Sang-Ho's notch filter script to remove power at frequencies corresponding to 1 year and 6 months.
    We are also removing a linear trend in this step.

    :param Data: a time-series object
    :returns: two time-series objects
    """

    Data = Data.remove_nans()

    # Parameters
    # %   x       1-D signal array
    # %   fs      sampling frequency, Hz
    # %   fn      notch frequency, Hz
    # %   Bn      notch bandwidth, Hz
    dt_interval = 1.0  # one day
    fs = 1 / dt_interval
    fn1 = 1.0 / 365.24  # fn = notch frequency, annual
    Bn1 = 0.1 * fn1
    fn2 = 2.0 / 365.24  # fn = notch frequency, semiannual
    Bn2 = 0.1 * fn2  # a choice: 10% seems to work well.

    # East, North, Up
    dE_filt = notch_filter.notchfilt(Data.dE, fs, fn1, Bn1, filtfiltopt=True)
    dE_filt = notch_filter.notchfilt(dE_filt, fs, fn2, Bn2, filtfiltopt=True)
    dE_detrended = simple_detrend_ts(Data.dtarray, dE_filt)

    dN_filt = notch_filter.notchfilt(Data.dN, fs, fn1, Bn1, filtfiltopt=True)
    dN_filt = notch_filter.notchfilt(dN_filt, fs, fn2, Bn2, filtfiltopt=True)
    dN_detrended = simple_detrend_ts(Data.dtarray, dN_filt)

    dU_filt = notch_filter.notchfilt(Data.dU, fs, fn1, Bn1, filtfiltopt=True)
    dU_filt = notch_filter.notchfilt(dU_filt, fs, fn2, Bn2, filtfiltopt=True)
    dU_detrended = simple_detrend_ts(Data.dtarray, dU_filt)

    detrended = Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=dN_detrended, dE=dE_detrended,
                           dU=dU_detrended, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes)
    trended = Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=dN_filt, dE=dE_filt, dU=dU_filt,
                         Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes)
    return detrended, trended


def remove_seasonals_by_STL(Data, STL_dir):
    """
    Has an issue: Right now only returns de-trended data.
    """
    filename = STL_dir + Data.name + "_STL_30.txt"

    if os.path.isfile(filename):
        # If a precomputed file exists...
        [dtstrings, dE, dN, dU, Se, Sn, Su] = np.loadtxt(filename, unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6),
                                                         dtype={'names': ('dt', 'dE', 'dN', 'dU', 'Se', 'Sn', 'Su'),
                                                                'formats': ('U8', np.float, np.float, np.float,
                                                                            np.float, np.float, np.float)})
        final_dtarray = [dt.datetime.strptime(x, "%Y%m%d") for x in dtstrings]
        Data = Timeseries(name=Data.name, coords=Data.coords, dtarray=final_dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se,
                          Su=Su, EQtimes=Data.EQtimes)

    else:  # ELSE: WE NEED TO RECOMPUTE
        print("Warning! STL not found for %s" % Data.name)
        print("We did not find a pre-computed array, so we are re-computing STL. ")

        # Preprocess data: remove nans, fill in gaps.
        Data = Data.remove_nans()
        [_, dE, Se] = preprocess_stl(Data.dtarray, Data.dE, Data.Se)
        [_, dN, Sn] = preprocess_stl(Data.dtarray, Data.dN, Data.Sn)
        [new_dtarray, dU, Su] = preprocess_stl(Data.dtarray, Data.dU, Data.Su)

        # Write E, N, U
        with open('raw_ts_data.txt', 'w') as ofile:
            for i in range(len(dE)):
                ofile.write('%s %f %f %f\n' % (dt.datetime.strftime(new_dtarray[i], "%Y%m%d"), dE[i], dN[i], dU[i]))

        # Call driver in matlab (read, STL, write)
        subprocess.call(['matlab', '-nodisplay', '-nosplash', '-r', 'stl_driver'], shell=False)

        # Read / Detrend the filtered data
        [dE, dN, dU] = np.loadtxt('filtered_ts_data.txt', unpack=True, usecols=(1, 2, 3))
        dE_detrended = simple_detrend_ts(new_dtarray, dE)
        dN_detrended = simple_detrend_ts(new_dtarray, dN)
        dU_detrended = simple_detrend_ts(new_dtarray, dU)

        # Put the gaps back in:
        final_dtarray, final_dE, final_dN, final_dU, final_Se, final_Sn, final_Su = [], [], [], [], [], [], []
        for i in range(len(new_dtarray)):
            if new_dtarray[i] in Data.dtarray:
                final_dtarray.append(new_dtarray[i])
                final_dE.append(dE_detrended[i])
                final_dN.append(dN_detrended[i])
                final_dU.append(dU_detrended[i])
                final_Se.append(Se[i])
                final_Sn.append(Sn[i])
                final_Su.append(Su[i])

        Data = Timeseries(name=Data.name, coords=Data.coords, dtarray=final_dtarray, dN=np.array(final_dN),
                          dE=np.array(final_dE), dU=np.array(final_dU), Sn=np.array(final_Sn), Se=np.array(final_Se),
                          Su=np.array(final_Su), EQtimes=Data.EQtimes)

        # Write the file so that we don't recompute it next time.
        io_other.write_stl(Data, filename)

    return Data, Data  # sorry, what is this?


def preprocess_stl(dtarray, data_column, uncertainties):
    """fill in gaps, and make to a multiple of 365 day cycles."""

    new_data_column, new_sig, new_dtarray = [], [], []
    new_data_column.append(data_column[0])
    new_sig.append(uncertainties[0])
    new_dtarray.append(dtarray[0])

    for i in range(1, len(dtarray)):

        start_counter = new_dtarray[-1]  # this is the date where we start.
        destination_counter = dtarray[i]  # the next datapoint we are going to append until in this loop.

        while start_counter + dt.timedelta(days=1) <= destination_counter:

            # If your next element is consecutive with the old element.
            if start_counter + dt.timedelta(days=1) == destination_counter:
                new_dtarray.append(dtarray[i])
                if data_column[i] == np.nan:
                    new_data_column.append(new_data_column[-1])
                    new_sig.append(uncertainties[-1])
                else:
                    new_data_column.append(data_column[i])
                    new_sig.append(uncertainties[i])
                break

            # If your next element is not consecutive with the old element
            else:
                new_dtarray.append(start_counter + dt.timedelta(days=1))
                new_data_column.append(new_data_column[-1])
                new_sig.append(uncertainties[-1])
                start_counter = start_counter + dt.timedelta(days=1)

    # Make the length an integer number of years
    while np.mod(len(new_dtarray), 365) > 0.01:
        new_dtarray.append(new_dtarray[-1] + dt.timedelta(days=1))
        new_data_column.append(new_data_column[-1])
        new_sig.append(uncertainties[-1])

    return [new_dtarray, new_data_column, new_sig]


def remove_seasonals_by_hydro(Data, hydro_dir, scaling=False):
    station = Data.name
    files = glob.glob(hydro_dir + station.lower() + '*.hyd')
    if not files:  # found an empty array
        print("Error! Hydro file not found for %s" % Data.name)
        print("Returning placeholder object.")
        wimpyObj = get_wimpy_object(Data)
        return wimpyObj, wimpyObj
    else:
        filename = files[0]

    # Read the hydro model, clean up the data, and pair it to the GPS
    hydro_data = io_nota.read_pbo_hydro_file(filename)
    Data = Data.remove_nans()
    hydro_data = hydro_data.remove_nans()
    [gps_data, hydro_data] = gps_ts_functions.pair_gps_model(Data, hydro_data)  # matched dtarrays.

    if scaling is True:
        [_, _, vert_gps] = gps_data.get_linear_annual_semiannual(critical_len=365)
        [_, _, vert_hydro] = hydro_data.get_linear_annual_semiannual(critical_len=365)
        gps_amp = np.sqrt(vert_gps[1] * vert_gps[1] + vert_gps[2] * vert_gps[2])
        hydro_amp = np.sqrt(vert_hydro[1] * vert_hydro[1] + vert_hydro[2] * vert_hydro[2])
        if hydro_amp == 0.0:
            print("ERROR! NLDAS amplitude is exactly 0!!  You should probably fix this. ")
            wimpyObj = get_wimpy_object(Data)
            print("Returning placeholder object.")
            return wimpyObj, wimpyObj
        scale_factor = gps_amp / hydro_amp
        print("NLDAS scaling factor is %.2f" % scale_factor)
    else:
        scale_factor = 1

    #  Subtract the model from the data.
    dE_filt = np.subtract(gps_data.dE, scale_factor*hydro_data.dE)
    dN_filt = np.subtract(gps_data.dN, scale_factor*hydro_data.dN)
    dU_filt = np.subtract(gps_data.dU, scale_factor*hydro_data.dU)

    # A Simple detrending
    dE_detrended = simple_detrend_ts(gps_data.dtarray, dE_filt)
    dN_detrended = simple_detrend_ts(gps_data.dtarray, dN_filt)
    dU_detrended = simple_detrend_ts(gps_data.dtarray, dU_filt)

    corrected_object = Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_detrended,
                                  dN=dN_detrended, dU=dU_detrended, Se=gps_data.Se,  Sn=gps_data.Sn, Su=gps_data.Su,
                                  EQtimes=gps_data.EQtimes)
    trended = Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_filt, dN=dN_filt,
                         dU=dU_filt, Se=gps_data.Se, Sn=gps_data.Sn, Su=gps_data.Su, EQtimes=gps_data.EQtimes)
    return corrected_object, trended


def remove_seasonals_by_german_load(Data, lsdm_dir):
    station = Data.name
    files = glob.glob(lsdm_dir + station + '*.txt')
    if not files:  # found an empty array
        print("Error! LSDM file not found for %s" % Data.name)
        print("Returning placeholder object.")
        wimpyObj = get_wimpy_object(Data)
        return wimpyObj, wimpyObj
    else:
        filename = files[0]

    # Read the hydro model, clean up data, and pair it to the GPS
    Data = Data.remove_nans()
    hydro_data = io_other.read_lsdm_file(filename)
    hydro_data = hydro_data.remove_nans()  # this may or may not be necessary.
    [gps_data, hydro_data] = gps_ts_functions.pair_gps_model(Data, hydro_data)  # matched in terms of dtarray.

    #  Subtract the model from the data.
    dE_filt = np.subtract(gps_data.dE, hydro_data.dE)
    dN_filt = np.subtract(gps_data.dN, hydro_data.dN)
    dU_filt = np.subtract(gps_data.dU, hydro_data.dU)

    # A simple detrending
    dE_detrended = simple_detrend_ts(gps_data.dtarray, dE_filt)
    dN_detrended = simple_detrend_ts(gps_data.dtarray, dN_filt)
    dU_detrended = simple_detrend_ts(gps_data.dtarray, dU_filt)

    detrended = Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_detrended,
                           dN=dN_detrended, dU=dU_detrended, Se=gps_data.Se, Sn=gps_data.Sn, Su=gps_data.Su,
                           EQtimes=gps_data.EQtimes)
    trended = Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=gps_data.dtarray, dE=dE_filt, dN=dN_filt,
                         dU=dU_filt, Se=gps_data.Se, Sn=gps_data.Sn, Su=gps_data.Su, EQtimes=gps_data.EQtimes)
    return detrended, trended


def remove_seasonals_by_lakes(Data, lakes_dir, lake_name):
    station = Data.name
    files = glob.glob(lakes_dir + station + "_" + lake_name + "*.txt")
    if not files:  # found an empty array
        print("Error! Lake %s file not found for %s" % (lake_name, Data.name))
        print("Returning placeholder object.")
        wimpyObj = get_wimpy_object(Data)
        return wimpyObj, wimpyObj
    else:
        filename = files[0]

    print("Reading loading TS from %s " % filename)
    loading_ts = io_other.read_lake_loading_ts(filename)
    [GPS_paired, loading_paired] = gps_ts_functions.pair_gps_model(Data, loading_ts)
    dE = np.subtract(GPS_paired.dE, loading_paired.dE)
    dN = np.subtract(GPS_paired.dN, loading_paired.dN)
    dU = np.subtract(GPS_paired.dU, loading_paired.dU)

    # A simple detrending
    dE_detrended = simple_detrend_ts(GPS_paired.dtarray, dE)
    dN_detrended = simple_detrend_ts(GPS_paired.dtarray, dN)
    dU_detrended = simple_detrend_ts(GPS_paired.dtarray, dU)

    corrected_object = Timeseries(name=Data.name, coords=Data.coords, dtarray=GPS_paired.dtarray, dE=dE_detrended,
                                  dN=dN_detrended, dU=dU_detrended, Se=GPS_paired.Se, Sn=GPS_paired.Sn,
                                  Su=GPS_paired.Su, EQtimes=Data.EQtimes)
    trended = Timeseries(name=Data.name, coords=Data.coords, dtarray=GPS_paired.dtarray, dE=dE, dN=dN, dU=dU,
                         Se=GPS_paired.Se, Sn=GPS_paired.Sn, Su=GPS_paired.Su, EQtimes=Data.EQtimes)
    return corrected_object, trended


def get_wimpy_object(Data):
    """Generate a timeseries object that has only lists of nans in its data fields."""
    placeholder = np.full_like(Data.dtarray, np.nan, dtype=np.double)
    wimpyObj = Timeseries(name=Data.name, coords=Data.coords, dtarray=Data.dtarray, dN=placeholder, dE=placeholder,
                          dU=placeholder, Sn=Data.Sn, Se=Data.Se, Su=Data.Su, EQtimes=Data.EQtimes)
    return wimpyObj


def remove_seasonals_by_GRACE(Data, grace_dir):
    """
    Here we use pre-computed GRACE load model time series to correct the GPS time series.
    We recognize that the horizontals will be bad, and that the resolution of GRACE is coarse.
    For these reasons, this is not an important part of the analysis.
    Read and interpolate GRACE loading model
    Subtract the GRACE model and remove it + overall trend from the GNSS time series
    """

    filename = grace_dir + "scaled_" + Data.name + "_PREM_model_ts.txt"
    if not os.path.isfile(filename):
        print("Error! GRACE not found for %s" % Data.name)
        print("Returning placeholder object.")
        wimpyObj = get_wimpy_object(Data)
        return wimpyObj, wimpyObj

    # If the station has been pre-computed with GRACE:
    Data = Data.remove_nans()
    grace_model = io_other.read_grace(filename)
    my_paired_ts = grace_ts_functions.pair_GPSGRACE(Data, grace_model)

    # Subtract the GRACE object
    dE_filt = np.subtract(my_paired_ts.east, my_paired_ts.u)
    dN_filt = np.subtract(my_paired_ts.north, my_paired_ts.v)
    dU_filt = np.subtract(my_paired_ts.vert, my_paired_ts.w)

    # A simple detrending
    dE_detrended = simple_detrend_ts(my_paired_ts.dtarray, dE_filt)
    dN_detrended = simple_detrend_ts(my_paired_ts.dtarray, dN_filt)
    dU_detrended = simple_detrend_ts(my_paired_ts.dtarray, dU_filt)

    detrended = Timeseries(name=Data.name, coords=Data.coords, dtarray=my_paired_ts.dtarray, dN=dN_detrended,
                           dE=dE_detrended, dU=dU_detrended, Sn=my_paired_ts.N_err, Se=my_paired_ts.E_err,
                           Su=my_paired_ts.V_err, EQtimes=Data.EQtimes)
    trended = Timeseries(name=Data.name, coords=Data.coords, dtarray=my_paired_ts.dtarray, dN=dN_filt, dE=dE_filt,
                         dU=dU_filt, Sn=my_paired_ts.N_err, Se=my_paired_ts.E_err, Su=my_paired_ts.V_err,
                         EQtimes=Data.EQtimes)
    return detrended, trended  # 0 = successful completion
