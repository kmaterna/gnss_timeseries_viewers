import numpy as np
from scipy.optimize import curve_fit


def linear_annual_semiannual_function(decyear, fit_params):
    """
    Build the function y = f(x), given curve parameters and set of observation times.
    Model consists of GPS_V = E*t + Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt)
    """
    model_def = []
    w = 2 * np.pi / 1.0
    for t in decyear:
        model_def.append(fit_params[0] * t + (fit_params[1] * np.cos(w * t)) + (fit_params[2] * np.sin(w * t)) + (
                fit_params[3] * np.cos(2 * w * t)) + (fit_params[4] * np.sin(2 * w * t)))
    return model_def


def annual_semiannual_only_function(decyear, fit_params):
    """
    Build the function y = f(x), given curve parameters and set of observation times.
    Model consists of GPS_V = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt)
    """
    model_def = []
    w = 2 * np.pi / 1.0
    for t in decyear:
        model_def.append(
            (fit_params[0] * np.cos(w * t)) + (fit_params[1] * np.sin(w * t)) + (fit_params[2] * np.cos(2 * w * t)) + (
                    fit_params[3] * np.sin(2 * w * t)))
    return model_def


def annual_only_function(decyear, fit_params):
    """
    Build the function y = f(x), given curve parameters and set of observation times.
    Model consists of GPS_V = Acos(wt) + Bsin(wt)
    """
    model_def = []
    w = 2 * np.pi / 1.0
    for t in decyear:
        model_def.append((fit_params[0] * np.cos(w * t)) + (fit_params[1] * np.sin(w * t)))
    return model_def


def construct_log_function(decday, fit_params):
    """
    Function has functional form:
    y = b + a*np.log(1+t/tau)
    fit params = [a, b, tau]
    The x axis is the same as the get_log_function()
    """
    a = fit_params[0]
    b = fit_params[1]
    tau = fit_params[2]
    model_def = []
    for i in range(len(decday)):
        model_def.append(b + a * np.log(1 + decday[i] / tau))
    return model_def


def invert_log_function(decyear, data_array):
    """
    Solve for function y = B + Alog(1+t/tau) relative to the day of the earthquake
    Same function as construct log function
    """
    def func(t, a, b, tau):
        return b + a * np.log(1 + t / tau)
    response = curve_fit(func, decyear, data_array)
    return response[0]


def invert_linear_annual_semiannual(decyear, data_array):
    """
    Fit a best-fitting equation to a time series in a linear least squares sense:
    GPS = Acos(wt) + Bsin(wt) + Ccos(2wt) + Dsin(2wt) + E*t + F
    Here we also solve for a linear trend as well.
    """
    design_matrix = []
    w = 2 * np.pi / 1.0
    for t in decyear:
        design_matrix.append([np.cos(w * t), np.sin(w * t), np.cos(2 * w * t), np.sin(2 * w * t), t, 1])
    design_matrix = np.array(design_matrix)
    params = np.dot(np.linalg.inv(np.dot(design_matrix.T, design_matrix)), np.dot(design_matrix.T, data_array))
    return params
