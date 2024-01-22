"""
Toolbox for operating on Timeseries objects.
Contains functions to map, filter, reduce, and process generic GPS time series
"""

import numpy as np
import datetime as dt
from scipy import signal
import scipy
import sys
from . import lssq_model_errors, utilities, math_functions
from Tectonic_Utils.geodesy import insar_vector_functions, haversine


class Timeseries:
    # The classic timeseries internal object!
    def __init__(self, name, coords, dtarray, dN, dE, dU, Sn, Se, Su, EQtimes=(), Offsettimes=()):
        self.name = name
        self.coords = coords
        self.dtarray = dtarray  # 1d array of datetime objects
        self.dN = dN  # 1d arrays in mm
        self.dE = dE
        self.dU = dU
        self.Sn = Sn  # 1d arrays in mm
        self.Se = Se
        self.Su = Su
        self.EQtimes = EQtimes  # just metadata, timing of earthquake offsets
        self.Offsettimes = Offsettimes  # just metadata, timing of antenna offsets
        if len(dtarray) == 0:  # zero-length timeseries throws error
            print("Error: length of dtarray is 0 for station %s. No timeseries created." % self.name)
            sys.exit(0)
        if len(dtarray) != len(dE) or len(dtarray) != len(dN) or len(dtarray) != len(dU):
            raise ValueError("Error! Length of data arrays and dt-array does not match!")
        if len(dtarray) != len(Se) or len(dtarray) != len(Sn) or len(dtarray) != len(Su):
            raise ValueError("Error! Length of uncertainty arrays and dt-array does not match!")

    # -------------------------------------------- #
    #              PREDICATES
    # -------------------------------------------- #

    def covers_date(self, include_time: dt.datetime) -> bool:
        if self.dtarray[0] < include_time < self.dtarray[-1]:
            return True
        else:
            return False

    def covers_date_range(self, include_time_start: dt.datetime, include_time_end: dt.datetime) -> bool:
        if self.dtarray[0] < include_time_start and self.dtarray[-1] > include_time_end:
            return True
        else:
            return False

    def has_incompatible_subwindow(self, starttime_desired: dt.datetime, endtime_desired: dt.datetime) -> bool:
        """
        Check for nasty things that disturb the calculation of slopes from sub-windowed timeseries,
        such as an incompatible window requested.
        """
        if endtime_desired < self.get_starttime():
            print("Error: end time before start of array for station %s. Returning Nan" % self.name)
            return True
        if starttime_desired > self.get_endtime():
            print("Error: start time after end of array for station %s. Returning Nan" % self.name)
            return True
        if starttime_desired >= endtime_desired:
            print("Error! Starttime after endtime for station %s " % self.name)
            return True
        if self.get_starttime() == self.get_endtime():
            print("Error: Timeseries has length 1")
            return True
        return False

    def is_within_bbox(self, bbox) -> bool:
        if bbox[0] <= self.coords[0] <= bbox[1] and bbox[2] <= self.coords[1] <= bbox[3]:
            return True
        else:
            return False

    # -------------------------------------------- #
    # PROCESSES THAT RETURN NEW TIME SERIES OBJECTS
    # -------------------------------------------- #

    def remove_outliers(self, outliers_def: float):
        medfilt_e = signal.medfilt(self.dE, 35)
        medfilt_n = signal.medfilt(self.dN, 35)
        medfilt_u = signal.medfilt(self.dU, 35)
        newdt, newdN, newdE, newdU, newSe, newSn, newSu = [], [], [], [], [], [], []
        for i in range(len(medfilt_e)):
            if abs(self.dE[i] - medfilt_e[i]) < outliers_def and abs(self.dN[i] - medfilt_n[i]) < outliers_def and abs(
                    self.dU[i] - medfilt_u[i]) < outliers_def * 2:
                newdt.append(self.dtarray[i])
                newdE.append(self.dE[i])
                newdN.append(self.dN[i])
                newdU.append(self.dU[i])
                newSe.append(self.Se[i])
                newSn.append(self.Sn[i])
                newSu.append(self.Su[i])
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=newdt, dN=np.array(newdN), dE=np.array(newdE),
                             dU=np.array(newdU), Sn=np.array(newSn), Se=np.array(newSe), Su=np.array(newSu),
                             EQtimes=self.EQtimes, Offsettimes=self.Offsettimes)
        return newData

    def impose_time_limits(self, starttime: dt.datetime, endtime: dt.datetime):
        """
        :param starttime: datetime object
        :param endtime:  datetime object
        :return: a Timeseries object
        """
        newdtarray, newdN, newdE, newdU, newSe, newSn, newSu = [], [], [], [], [], [], []
        for i in range(len(self.dN)):
            if starttime <= self.dtarray[i] <= endtime:
                newdtarray.append(self.dtarray[i])
                newdE.append(self.dE[i])
                newdN.append(self.dN[i])
                newdU.append(self.dU[i])
                newSe.append(self.Se[i])
                newSn.append(self.Sn[i])
                newSu.append(self.Su[i])
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=newdtarray, dN=np.array(newdN),
                             dE=np.array(newdE), dU=np.array(newdU), Sn=np.array(newSn), Se=np.array(newSe),
                             Su=np.array(newSu), EQtimes=self.EQtimes, Offsettimes=self.Offsettimes)
        return newData

    def remove_nans(self):
        idxE, idxN, idxU = np.isnan(self.dE), np.isnan(self.dN), np.isnan(self.dU)
        temp_dates, temp_east, temp_north, temp_vert = [], [], [], []
        temp_Sn, temp_Se, temp_Su = [], [], []
        if (sum(idxN) + sum(idxE) + sum(idxU)) == 0:
            return self
        else:  # if there are nans, please pull them out.
            for i in range(len(self.dtarray)):
                if idxE[i] == 0 and idxN[i] == 0 and idxU[i] == 0:
                    temp_dates.append(self.dtarray[i])
                    temp_east.append(self.dE[i])
                    temp_north.append(self.dN[i])
                    temp_vert.append(self.dU[i])
                    temp_Se.append(self.Se[i])
                    temp_Sn.append(self.Sn[i])
                    temp_Su.append(self.Su[i])
            newData = Timeseries(name=self.name, coords=self.coords, dtarray=temp_dates, dN=np.array(temp_north),
                                 dE=np.array(temp_east), dU=np.array(temp_vert), Sn=np.array(temp_Sn),
                                 Se=np.array(temp_Se), Su=np.array(temp_Su), EQtimes=self.EQtimes,
                                 Offsettimes=self.Offsettimes)
        return newData

    def remove_specific_date(self, target_date: dt.datetime):
        """
        Exclude a particular date from a time series. Useful when needing to exclude an earthquake event day.

        :param target_date: datetime object
        :return: Timeseries
        """
        newdtarray, newdN, newdE, newdU, newSe, newSn, newSu = [], [], [], [], [], [], []
        for i in range(len(self.dN)):
            if self.dtarray[i] != target_date:
                newdtarray.append(self.dtarray[i])
                newdE.append(self.dE[i])
                newdN.append(self.dN[i])
                newdU.append(self.dU[i])
                newSe.append(self.Se[i])
                newSn.append(self.Sn[i])
                newSu.append(self.Su[i])
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=newdtarray, dN=np.array(newdN),
                             dE=np.array(newdE), dU=np.array(newdU), Sn=np.array(newSn), Se=np.array(newSe),
                             Su=np.array(newSu), EQtimes=self.EQtimes, Offsettimes=self.Offsettimes)
        return newData

    def remove_constant(self, east_offset=0.0, north_offset=0.0, vert_offset=0.0):
        """Subtract a constant number from each data array in a time series object

        :param east_offset: a scalar, in mm
        :param north_offset: a scalar, in mm
        :param vert_offset: a scalar, in mm
        """
        temp_east = np.array([x - east_offset for x in self.dE])
        temp_north = np.array([x - north_offset for x in self.dN])
        temp_vert = np.array([x - vert_offset for x in self.dU])
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=self.dtarray,
                             dN=np.array(temp_north), dE=np.array(temp_east), dU=np.array(temp_vert),
                             Sn=self.Sn, Se=self.Se, Su=self.Su, EQtimes=self.EQtimes, Offsettimes=self.Offsettimes)
        return newData

    def set_zero_mean(self):
        """
        Subtract a constant such that each of dE, dN, and dU have zero mean.
        """
        east0 = float(np.nanmean(self.dE))
        north0 = float(np.nanmean(self.dN))
        vert0 = float(np.nanmean(self.dU))
        newData = self.remove_constant(east_offset=east0, north_offset=north0, vert_offset=vert0)
        return newData

    def embed_tsobject_with_eqdates(self, eq_obj):
        """
        :param eq_obj: a list of Offset objects to embed into the Timeseries earthquake metadata
        """
        eqdates = [x.evdt for x in eq_obj]
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=self.dtarray, dN=self.dN,
                             dE=self.dE, dU=self.dU, Sn=self.Sn, Se=self.Se,
                             Su=self.Su, EQtimes=eqdates, Offsettimes=self.Offsettimes)
        return newData

    def embed_tsobject_with_offset_dates(self, offset_obj):
        """
        :param offset_obj: a list of Offset objects to embed into the Timeseries offset metadata
        """
        offsetdates = [x.evdt for x in offset_obj]
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=self.dtarray, dN=self.dN,
                             dE=self.dE, dU=self.dU, Sn=self.Sn, Se=self.Se,
                             Su=self.Su, EQtimes=self.EQtimes, Offsettimes=offsetdates)
        return newData

    def rotate_data(self, azimuth: float):  # should test to make sure this works the way I think it does.
        """
        :param azimuth: degrees, CW from north
        :return: a Timeseries object. ts.e is associated with the new azimuth ts.e is associated with azimuth+90.
        """
        newE, newN = [], []
        theta = insar_vector_functions.bearing_to_cartesian(azimuth)
        for i in range(len(self.dtarray)):
            new_position = insar_vector_functions.rotate_vector_by_angle(self.dE[i], self.dN[i], theta)
            newE.append(new_position[0])
            newN.append(new_position[1])
        rotated_ts = Timeseries(name=self.name, coords=self.coords, dtarray=self.dtarray,
                                dN=np.array(newN), dE=np.array(newE), dU=self.dU,
                                Sn=self.Sn, Se=self.Se, Su=self.Su, EQtimes=self.EQtimes, Offsettimes=self.Offsettimes)
        return rotated_ts

    def median_filter_data(self, size: int):
        """
        :param size: integer, window size for median filter
        :return: a Timeseries object
        """
        udata = scipy.ndimage.median_filter(self.dE, size=size)
        vdata = scipy.ndimage.median_filter(self.dN, size=size)
        wdata = scipy.ndimage.median_filter(self.dU, size=size)
        filtered_ts = Timeseries(dtarray=self.dtarray, dE=udata, dN=vdata, dU=wdata, coords=self.coords,
                                 name=self.name, Se=self.Se, Sn=self.Sn, Su=self.Su, EQtimes=self.EQtimes,
                                 Offsettimes=self.Offsettimes)
        return filtered_ts

    def detrend_data_by_value(self, east_params, north_params, vert_params):
        """
        :param east_params: list of floats [slope, a2(cos), a1(sin), s2(cos) s2(sin)]
        :param north_params: list of floats [slope, a2(cos), a1(sin), s2(cos) s2(sin)]
        :param vert_params: list of floats [slope, a2(cos), a1(sin), s2(cos) s2(sin)]
        :return: a detrended Timeseries object
        """
        if sum(np.isnan(east_params)) > 0 or sum(np.isnan(north_params)) > 0 or sum(np.isnan(vert_params)) > 0:
            print("ERROR: Your input slope values contain nan!")
            return self

        # Parameters Format: slope, a2(cos), a1(sin), s2, s1.
        east_detrended, north_detrended, vert_detrended = [], [], []
        data = self.remove_nans()
        decyear = utilities.get_float_times(data.dtarray)

        east_model = math_functions.linear_annual_semiannual_function(decyear, east_params)
        north_model = math_functions.linear_annual_semiannual_function(decyear, north_params)
        vert_model = math_functions.linear_annual_semiannual_function(decyear, vert_params)

        for i in range(len(decyear)):
            east_detrended.append(data.dE[i] - (east_model[i]))
            north_detrended.append(data.dN[i] - (north_model[i]))
            vert_detrended.append(data.dU[i] - (vert_model[i]))
        east_detrended = [x - east_detrended[0] for x in east_detrended]
        north_detrended = [x - north_detrended[0] for x in north_detrended]
        vert_detrended = [x - vert_detrended[0] for x in vert_detrended]
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=data.dtarray, dN=north_detrended,
                             dE=east_detrended, dU=vert_detrended, Sn=data.Sn, Se=data.Se, Su=data.Su,
                             EQtimes=data.EQtimes, Offsettimes=self.Offsettimes)
        return newData

    def simple_detrend(self, starttime=None, endtime=None):
        """
        The simplest de-trending option just based upon the linear slope from starttime to endttime.

        :param starttime: dt.datetime
        :param endtime: dt.datetime
        :return: a Timeseries object
        """
        [e_slope, n_slope, u_slope, _, _, _] = self.get_slope(starttime=starttime, endtime=endtime)
        newData = self.detrend_data_by_value([e_slope, 0, 0, 0, 0], [n_slope, 0, 0, 0, 0], [u_slope, 0, 0, 0, 0])
        return newData

    def remove_seasonal_by_value(self, east_params, north_params, vert_params):
        """
        Least squares seasonal parameters. Remove seasonal components.
        Parameters Format: slope, a2(cos), a1(sin), s2, s1.
        """
        east_detrended, north_detrended, vert_detrended = [], [], []
        data = self.remove_nans()
        decyear = utilities.get_float_times(data.dtarray)

        east_model = math_functions.annual_semiannual_only_function(decyear, east_params[1:])
        north_model = math_functions.annual_semiannual_only_function(decyear, north_params[1:])
        vert_model = math_functions.annual_semiannual_only_function(decyear, vert_params[1:])

        for i in range(len(decyear)):
            east_detrended.append(data.dE[i] - (east_model[i]))
            north_detrended.append(data.dN[i] - (north_model[i]))
            vert_detrended.append(data.dU[i] - (vert_model[i]))
        newData = Timeseries(name=self.name, coords=self.coords, dtarray=data.dtarray, dN=north_detrended,
                             dE=east_detrended, dU=vert_detrended, Sn=data.Sn, Se=data.Se, Su=data.Su,
                             EQtimes=data.EQtimes, Offsettimes=self.Offsettimes)
        return newData

    # -------------------------------------------- #
    # REDUCE TIME SERIES OBJECTS TO SCALARS OR VALUES
    # -------------------------------------------- #

    def get_starttime(self) -> dt.datetime:
        return self.dtarray[0]

    def get_endtime(self) -> dt.datetime:
        return self.dtarray[-1]

    def get_time_range(self) -> (dt.datetime, dt.datetime):
        return self.get_starttime(), self.get_endtime()

    def get_distance_to_point(self, point) -> float:
        """
        :param point: (lon, lat)
        :return: distance in km
        """
        return haversine.distance((self.coords[1], self.coords[0]), (point[1], point[0]))

    def get_slope(self, starttime=None, endtime=None, missing_fraction=0.6):
        """
        Model data with a best-fit y = mx + b.
        Returns six numbers: e_slope, n_slope, v_slope, e_std, n_std, v_std

        :param starttime: dt datetime
        :param endtime: dt datetime
        :param missing_fraction: float between 0 and 1
        """
        if starttime is None:
            starttime = self.get_starttime()
        if endtime is None:
            endtime = self.get_endtime()

        # Defensive programming
        if self.has_incompatible_subwindow(starttime_desired=starttime, endtime_desired=endtime):
            print("Basic defensive programming failed. Returning nans for slopes")
            return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

        # Cut to desired window, and remove nans
        data = self.remove_nans()
        data = data.impose_time_limits(starttime, endtime)

        # More defensive programming
        if len(data.dtarray) <= 2:
            print("ERROR: no time array for station %s. Returning Nan" % data.name)
            return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        time_duration = data.dtarray[-1] - data.dtarray[0]
        if time_duration.days < 270:
            print(
                "ERROR: using <<<1 year of data to estimate parameters for station %s. Returning Nan" % data.name)
            return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        if len(data.dE) < time_duration.days * missing_fraction:
            print(
                "ERROR: Most of the data is missing to estimate parameters for station %s. Returning Nan" % data.name)
            return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

        # doing the inversion here, since it's only one line.
        decyear = utilities.get_float_times(data.dtarray)
        east_coef = np.polyfit(decyear, data.dE, 1)
        north_coef = np.polyfit(decyear, data.dN, 1)
        vert_coef = np.polyfit(decyear, data.dU, 1)
        east_slope, north_slope, vert_slope = east_coef[0], north_coef[0], vert_coef[0]

        # How bad is the fit to the line?
        east_trend = [east_coef[0] * x + east_coef[1] for x in decyear]
        east_detrended = [data.dE[i] - east_trend[i] for i in range(len(data.dE))]
        east_std = np.std(east_detrended)
        north_trend = [north_coef[0] * x + north_coef[1] for x in decyear]
        north_detrended = [data.dN[i] - north_trend[i] for i in range(len(data.dN))]
        north_std = np.std(north_detrended)
        vert_trend = [vert_coef[0] * x + vert_coef[1] for x in decyear]
        vert_detrended = [data.dU[i] - vert_trend[i] for i in range(len(data.dU))]
        vert_std = np.std(vert_detrended)

        return [east_slope, north_slope, vert_slope, east_std, north_std, vert_std]

    def get_slope_unc(self, starttime: dt.datetime, endtime: dt.datetime):
        """
        Calculate Allan Variance of Rates.
        slope is returned as params[0]
        """
        data = self.impose_time_limits(starttime, endtime)
        x = utilities.get_float_times(data.dtarray)
        params, covm = lssq_model_errors.AVR(x, data.dE, data.Se, verbose=0)
        Esigma = np.sqrt(covm[0][0])
        params, covm = lssq_model_errors.AVR(x, data.dN, data.Sn, verbose=0)
        Nsigma = np.sqrt(covm[0][0])
        params, covm = lssq_model_errors.AVR(x, data.dU, data.Su, verbose=0)
        Usigma = np.sqrt(covm[0][0])
        return [Esigma, Nsigma, Usigma]

    def get_mean_position(self, starttime=None, endtime=None):
        """
        Compute average value of the time series between starttime and endtime

        :param starttime: dt.datetime object, default is the start of the timeseries
        :param endtime: dt.datetime object, default is the end of the timeseries
        :returns: [E, N, U] list of floats
        """
        # Defensive programming
        if starttime is None:
            starttime = self.get_starttime()
        if endtime is None:
            endtime = self.get_endtime()

        if self.has_incompatible_subwindow(starttime_desired=starttime, endtime_desired=endtime):
            return [np.nan, np.nan, np.nan]
        data = self.remove_nans()  # Cut to desired window, and remove nans
        data = data.impose_time_limits(starttime, endtime)
        return [np.nanmean(data.dE), np.nanmean(data.dN), np.nanmean(data.dU)]

    def get_values_at_date(self, selected_date: dt.datetime, num_days=10):
        """
        At a selected date, extract east, north, and up values at that date, with a forward-looking averaging window.
        """
        if selected_date in self.dtarray:
            [e_value, n_value, u_value] = self.get_mean_position(starttime=selected_date,
                                                                 endtime=selected_date+dt.timedelta(days=num_days))
        else:
            print("Error: requested date %s not found in dtarray" % dt.datetime.strftime(selected_date, "%Y-%m-%d"))
            [e_value, n_value, u_value] = [np.nan, np.nan, np.nan]
        return e_value, n_value, u_value

    def subsample_in_time(self, target_time: dt.datetime, window_days=30):
        """
        Downsample the position corresponding a given date by averaging over a two-sided window around target date.
        Almost the same as the get_values_at_date(), just with two-sided averaging window.

        :param target_time: dt.datetime object
        :param window_days: int, length of two-sided window in days
        :returns: E0, N0, U0
        """
        dE_start, dN_start, dU_start = [], [], []
        for i in range(len(self.dtarray)):
            if abs((self.dtarray[i] - target_time).days) < window_days:
                dE_start.append(self.dE[i])
                dN_start.append(self.dN[i])
                dU_start.append(self.dU[i])
        if len(dE_start) > 2:
            E0, N0, U0 = np.nanmean(dE_start), np.nanmean(dN_start), np.nanmean(dU_start)
        else:
            E0, N0, U0 = np.nan, np.nan, np.nan
        return E0, N0, U0

    def get_uncertainties_for_avg_position(self, selected_date: dt.datetime, num_days: int):
        """
        Uncertainties on the average displacements are computed multiplying by 2/(sqrt(n)) because
        we approximate it as ~2x larger than uncorrelated noise.
        GPS uncertainties are underestimated by factors of 2-11 if time-correlated noise is not considered.
        Hackl [2011] Johnson and Agnew [1995] Zhang et al. [1997] Mao et al. [1999] Williams et al. [2004].
        Still a work in progress.
        """
        Se, Sn, Su = np.nan, np.nan, np.nan
        if selected_date in self.dtarray:
            idx = self.dtarray.index(selected_date)
            Se = self.Se[idx] * (2 / np.sqrt(num_days))
            Sn = self.Sn[idx] * (2 / np.sqrt(num_days))
            Su = self.Su[idx] * (2 / np.sqrt(num_days))
        else:
            print("Error: requested date %s not found in dtarray" % dt.datetime.strftime(selected_date, "%Y-%m-%d"))
        return Se, Sn, Su

    def get_gap_fraction(self) -> float:
        """Compute the fraction of the time series that is not present. """
        starttime = self.dtarray[0]
        endtime = self.dtarray[-1]
        delta = endtime - starttime
        total_possible_days = delta.days
        days_present = len(self.dtarray)
        return 1 - (days_present / total_possible_days)

    def get_logfunction_params(self, eqtime: dt.datetime):
        """
        y = B + Alog(1+t/tau), t is in decyear
        Useful for postseismic transients.

        :param eqtime: dt.datetime object
        :returns: [list of e_params, list of n_params, list of u_params]
        """
        float_times = utilities.get_relative_times(self.dtarray, eqtime)  # in days
        e_params, ecov = math_functions.invert_log_function(float_times, self.dE)
        n_params, ecov = math_functions.invert_log_function(float_times, self.dN)
        u_params, ecov = math_functions.invert_log_function(float_times, self.dU)
        return [e_params, n_params, u_params]

    def get_linear_annual_semiannual(self, starttime=None, endtime=None, critical_len=365):
        """
        The critical_len parameter allows us to switch this function for both GPS and GRACE time series in GPS format
        Model the data with a best-fit GPS = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt) + E*t + F

        :param starttime: dt.datetime object, default start of time series
        :param endtime: dt.datetime object, default end of time series
        :param critical_len: int, minimum duration of time series in days
        :returns: [list of e_params, list of n_params, list of u_params]
        """
        if starttime is None:
            starttime = self.get_starttime()
        if endtime is None:
            endtime = self.get_endtime()

        # Defensive programming
        if self.has_incompatible_subwindow(starttime_desired=starttime, endtime_desired=endtime):
            east_params = [np.nan, 0, 0, 0, 0]
            north_params = [np.nan, 0, 0, 0, 0]
            vert_params = [np.nan, 0, 0, 0, 0]
            return [east_params, north_params, vert_params]

        # Cut to desired window, and remove nans
        x = self.remove_nans()
        x = x.impose_time_limits(starttime, endtime)

        duration = x.dtarray[-1] - x.dtarray[0]
        if duration.days < critical_len:
            print(
                "ERROR: using less than 1 year of data to estimate params for station %s. Returning Nan" % x.name)
            east_params = [np.nan, 0, 0, 0, 0]
            north_params = [np.nan, 0, 0, 0, 0]
            up_params = [np.nan, 0, 0, 0, 0]
            return [east_params, north_params, up_params]

        decyear = utilities.get_float_times(x.dtarray)
        east_params_unordered = math_functions.invert_linear_annual_semiannual(decyear, x.dE)
        north_params_unordered = math_functions.invert_linear_annual_semiannual(decyear, x.dN)
        vert_params_unordered = math_functions.invert_linear_annual_semiannual(decyear, x.dU)

        # The definition for returning parameters:
        # slope, a2(cos), a1(sin), s2, s1.
        east_params = [east_params_unordered[4], east_params_unordered[0], east_params_unordered[1],
                       east_params_unordered[2], east_params_unordered[3]]
        north_params = [north_params_unordered[4], north_params_unordered[0], north_params_unordered[1],
                        north_params_unordered[2], north_params_unordered[3]]
        vert_params = [vert_params_unordered[4], vert_params_unordered[0], vert_params_unordered[1],
                       vert_params_unordered[2], vert_params_unordered[3]]

        return [east_params, north_params, vert_params]


# -------------------------------------------- #
# FUNCTIONS THAT OPERATE ON TWO TIMESERIES OBJECTS
# -------------------------------------------- #

def pair_gps_model(gps_data: Timeseries, model_data: Timeseries):
    """
    Take two time series objects, and return two paired time series objects.
    It could be that GPS has days that model doesn't, or the other way around.
    """
    dtarray, dE_gps, dN_gps, dU_gps, Se_gps, Sn_gps, Su_gps = [], [], [], [], [], [], []
    dE_model, dN_model, dU_model, Se_model, Sn_model, Su_model = [], [], [], [], [], []
    gps_data = gps_data.remove_nans()
    model_data = model_data.remove_nans()
    for i in range(len(gps_data.dtarray)):
        if gps_data.dtarray[i] in model_data.dtarray:
            idx = model_data.dtarray.index(gps_data.dtarray[i])  # where is this datetime object in the model array?
            dtarray.append(gps_data.dtarray[i])
            dE_gps.append(gps_data.dE[i])
            dN_gps.append(gps_data.dN[i])
            dU_gps.append(gps_data.dU[i])
            Se_gps.append(gps_data.Se[i])
            Sn_gps.append(gps_data.Sn[i])
            Su_gps.append(gps_data.Su[i])
            dE_model.append(model_data.dE[idx])
            dN_model.append(model_data.dN[idx])
            dU_model.append(model_data.dU[idx])
            Se_model.append(model_data.Se[idx])
            Sn_model.append(model_data.Sn[idx])
            Su_model.append(model_data.Su[idx])
    paired_gps = Timeseries(name=gps_data.name, coords=gps_data.coords, dtarray=dtarray, dE=dE_gps, dN=dN_gps,
                            dU=dU_gps, Se=Se_gps, Sn=Sn_gps, Su=Su_gps, EQtimes=gps_data.EQtimes,
                            Offsettimes=gps_data.Offsettimes)
    paired_model = Timeseries(name=model_data.name, coords=model_data.coords, dtarray=dtarray, dE=dE_model, dN=dN_model,
                              dU=dU_model, Se=Se_model, Sn=Sn_model, Su=Su_model, EQtimes=model_data.EQtimes,
                              Offsettimes=model_data.Offsettimes)
    return [paired_gps, paired_model]


def pair_gps_model_keeping_gps(gps_data: Timeseries, model_data: Timeseries):
    """
    Take two time series objects, and return two time series objects.
    Keep all data from the first one.
    Generate a model_data with length that matches the GPS.
    """
    dtarray, dE_model, dN_model, dU_model, Se_model, Sn_model, Su_model = [], [], [], [], [], [], []
    gps_data = gps_data.remove_nans()
    model_data = model_data.remove_nans()
    for i in range(len(gps_data.dtarray)):
        if gps_data.dtarray[i] in model_data.dtarray:
            idx = model_data.dtarray.index(gps_data.dtarray[i])  # where is this datetime object in the model array?
            dtarray.append(gps_data.dtarray[i])
            dE_model.append(model_data.dE[idx])
            dN_model.append(model_data.dN[idx])
            dU_model.append(model_data.dU[idx])
            Se_model.append(model_data.Se[idx])
            Sn_model.append(model_data.Sn[idx])
            Su_model.append(model_data.Su[idx])
        else:
            dtarray.append(gps_data.dtarray[i])  # if we can't find it, then we put filler model.
            dE_model.append(0)
            dN_model.append(0)
            dU_model.append(0)
            Se_model.append(0)
            Sn_model.append(0)
            Su_model.append(0)
    paired_model = Timeseries(name=model_data.name, coords=model_data.coords, dtarray=dtarray, dE=dE_model, dN=dN_model,
                              dU=dU_model, Se=Se_model, Sn=Sn_model, Su=Su_model, EQtimes=model_data.EQtimes,
                              Offsettimes=model_data.Offsettimes)
    return [gps_data, paired_model]


def get_referenced_data(roving_station_data: Timeseries, base_station_data: Timeseries) -> Timeseries:
    """
    Take a time series object and remove motion of a base station (another time series object)
    If there's a starttime, then we will solve for a best-fitting model offset at starttime.
    Used when subtracting models
    """
    dtarray, dE_gps, dN_gps, dU_gps, Se_gps, Sn_gps, Su_gps = [], [], [], [], [], [], []
    roving_station_data = roving_station_data.remove_nans()
    for i in range(len(roving_station_data.dtarray)):
        if roving_station_data.dtarray[i] in base_station_data.dtarray:
            idx = base_station_data.dtarray.index(
                roving_station_data.dtarray[i])  # where is this datetime object in the model array?
            dtarray.append(roving_station_data.dtarray[i])
            dE_gps.append(roving_station_data.dE[i] - base_station_data.dE[idx])
            dN_gps.append(roving_station_data.dN[i] - base_station_data.dN[idx])
            dU_gps.append(roving_station_data.dU[i] - base_station_data.dU[idx])
            Se_gps.append(roving_station_data.Se[i])
            Sn_gps.append(roving_station_data.Sn[i])
            Su_gps.append(roving_station_data.Su[i])
    gps_relative = Timeseries(name=roving_station_data.name, coords=roving_station_data.coords, dtarray=dtarray,
                              dE=np.array(dE_gps), dN=np.array(dN_gps), dU=np.array(dU_gps), Se=np.array(Se_gps),
                              Sn=np.array(Sn_gps), Su=np.array(Su_gps), EQtimes=roving_station_data.EQtimes,
                              Offsettimes=roving_station_data.Offsettimes)
    return gps_relative

# -------------------------------------------- #
# FOR REVERSE COMPATIBILITY (MAY BE DEPRECATED)
# -------------------------------------------- #

def remove_outliers(Data0: Timeseries, outliers_def: float):
    return Data0.remove_outliers(outliers_def)

def impose_time_limits(Data0: Timeseries, starttime: dt.datetime, endtime: dt.datetime):
    return Data0.impose_time_limits(starttime, endtime)

def remove_nans(Data0: Timeseries):
    return Data0.remove_nans()

def remove_constant(Data0: Timeseries, east_offset=0, north_offset=0, vert_offset=0):
    return Data0.remove_constant(east_offset, north_offset, vert_offset)

def detrend_data_by_value(Data0: Timeseries, east_params, north_params, vert_params):
    return Data0.detrend_data_by_value(east_params, north_params, vert_params)

def remove_seasonal_by_value(Data0: Timeseries, east_params, north_params, vert_params):
    return Data0.remove_seasonal_by_value(east_params, north_params, vert_params)

def get_slope(Data0: Timeseries, starttime=None, endtime=None, missing_fraction=0.6):
    return Data0.get_slope(starttime, endtime, missing_fraction)

def get_slope_unc(dataObj: Timeseries, starttime, endtime):
    return dataObj.get_slope_unc(starttime, endtime)

def get_linear_annual_semiannual(Data0: Timeseries, starttime=None, endtime=None, critical_len=365):
    return Data0.get_linear_annual_semiannual(starttime, endtime, critical_len)

def get_mean_position(Data0: Timeseries, starttime=None, endtime=None):
    return Data0.get_mean_position(starttime, endtime)

# -------------------------------------------- #
# FUNCTIONS ON LISTS OF TS OBJECTS
# -------------------------------------------- #

def get_starttime(ts_objects_list) -> dt.datetime:
    """
    :param ts_objects_list: a list of ts objects
    :return: earliest time of any ts object
    """
    starttime = dt.datetime.strptime('2050-01-01', "%Y-%m-%d")  # crazy guess
    for item in ts_objects_list:
        if item.dtarray[0] < starttime:
            starttime = item.dtarray[0]
    return starttime


def get_endtime(ts_objects_list) -> dt.datetime:
    """
    :param ts_objects_list: a list of ts objects
    :return: latest time of any ts object
    """
    endtime = dt.datetime.strptime('1950-01-01', "%Y-%m-%d")  # crazy guess
    for item in ts_objects_list:
        if item.dtarray[-1] > endtime:
            endtime = item.dtarray[-1]
    return endtime


def get_starttime_endtime(ts_objects_list):
    return get_starttime(ts_objects_list), get_endtime(ts_objects_list)


def filter_to_station_subset(ts_objects_list, desired_station_list):
    """
    :param ts_objects_list: list of Timeseries objects
    :param desired_station_list: list of strings
    :return: list of Timeseries objects
    """
    return [x for x in ts_objects_list if x.name in desired_station_list]


def get_mean_ts(list_of_ts_objects) -> Timeseries:
    """
    Take a number of timeseries objects and create a new one that's the mean of the originals.
    """
    list_of_ts_objects = [x.set_zero_mean() for x in list_of_ts_objects]
    starttime = get_starttime(list_of_ts_objects)
    endtime = get_endtime(list_of_ts_objects)
    dtarray, dE_mean, dN_mean, dU_mean, Se_mean, Sn_mean, Su_mean = [], [], [], [], [], [], []
    target_date = starttime
    while target_date < endtime:
        dtarray.append(target_date)
        east_vals, north_vals, up_vals = [], [], []
        for station in list_of_ts_objects:
            if target_date in station.dtarray:
                e_value, n_value, u_value = station.get_values_at_date(target_date, num_days=1)
                east_vals.append(e_value)
                north_vals.append(n_value)
                up_vals.append(u_value)
        dE_mean.append(np.nanmean(east_vals))
        dN_mean.append(np.nanmean(north_vals))
        dU_mean.append(np.nanmean(up_vals))
        Se_mean.append(np.nanstd(east_vals))
        Sn_mean.append(np.nanstd(north_vals))
        Su_mean.append(np.nanstd(up_vals))
        target_date = target_date + dt.timedelta(days=1)
    avg_ts = Timeseries(name='mean', coords=(np.nan, np.nan), dtarray=dtarray,
                        dE=np.array(dE_mean), dN=np.array(dN_mean), dU=np.array(dU_mean), Se=np.array(Se_mean),
                        Sn=np.array(Sn_mean), Su=np.array(Su_mean), EQtimes=list_of_ts_objects[0].EQtimes,
                        Offsettimes=list_of_ts_objects[0].Offsettimes)
    return avg_ts
