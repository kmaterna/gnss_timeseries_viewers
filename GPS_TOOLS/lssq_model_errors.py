# A toolbox for fitting by least squares and getting uncertainties on velocity parameters. 
# TOOLS:
# Linear_fitting_menke  : the formal uncertainties from the textbook
# fit_curvefit          : a python function that does practically the same thing, but with more sophisticated handling of individual errors
# colored_noise         : Allan Variance

import numpy as np 
from scipy import optimize


def linear_fitting_menke(x, y, sig, verbose=1):
	"""
	Take a time series and fit a best-fitting linear least squares equation: 
	GPS = M*t + B; 
	We get the sigma on the model parameters from Page 64 of Menke, chapter 3
	cov(m) = sigma^2 * [G^T G]^-1
	In this case, sigma is a uniform value for each data point
	
	I don't really like this option, since it assumes random white noise and only cares about the 
	uncertainties and x-placement of the data (doesn't care what the data actually is). 
	However, it's easy, and technically correct. 
	"""
	if verbose:
		print("\nEstimating slope and intercept from Menke Error Propagation");
	design_matrix=[];
	for t in x:
		design_matrix.append([t, 1]);
	design_matrix= np.array(design_matrix);
	params = np.dot(np.linalg.inv(np.dot(design_matrix.T, design_matrix)), np.dot(design_matrix.T, y));
	covm = np.linalg.inv(np.dot(design_matrix.T, design_matrix));
	covm = np.multiply(covm, sig*sig);

	error = [] 
	for i in range(len(params)):
		try:
			error.append(np.absolute(covm[i][i])**0.5)
		except:
			error.append( 0.00 )
	perr = np.array(error)
	if verbose:
		print("pfit = ", params);
		print("perr = ", perr)

	return params, covm;


def fit_curvefit(x, y, sig, verbose=1):
	"""
	Note: As per the current documentation (Scipy V1.1.0), sigma (yerr) must be:
	None or M-length sequence or MxM array, optional
	p0 is an initial guess for the first parameter
	f is the function we're fitting (line in this case)
	sigma must be an array
	"""

	if verbose:
		print("\nEstimating slope and intercept from scipy.optimize.curve_fit method");
	def f(x, p0, p1):
		return p0*x + p1;
		# The function we are curve fitting

	pfit, pcov = optimize.curve_fit(f,x,y, sigma=sig, absolute_sigma=True, epsfcn=0.0001);
	error = [] 
	for i in range(len(pfit)):
		try:
			error.append(np.absolute(pcov[i][i])**0.5)
		except:
			error.append( 0.00 )
	perr = np.array(error)
	
	if verbose:
		print("# Fit parameters and parameter errors from curve_fit method :")
		print("pfit = ", pfit);
		print("perr = ", perr);

	return pfit, pcov;


def AVR(x, y, sig, verbose=1, overlapping=True):
	# This is based on the Allan Variance for Rates, devised in Hackl et al. 
	# (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JB008142)
	# This has 2 tunable parameters, tau and the overlapping inveral
	# This provides very similar uncertainty estimates to MLE techniques, but it's faster. 

	tau = 150;
	step_interval = 3;  # how much do you step to have overlapping windows
	
	if overlapping==False:
		step_interval=tau;  # non-overlapping windows case

	if tau > len(x)*0.25:
		print("Error! Tau was more than 1/4 of the dataset. Trimming tau to len(x)*0.25. ");
		tau=int(np.floor(len(x)*0.25));

	if len(x)<10:
		print("Error! too short time series. No AVR computed");
		return [np.nan, np.nan], [[np.nan, np.nan], [np.nan, np.nan]];
	
	if verbose:
		print("\nEstimating slope and intercept from AVR method with tau = %d" % tau);		

	slopes=[];
	params_overall, cov_overall = linear_fitting_menke(x, y, np.mean(sig), verbose=0);
	for i in range(0,len(x)-tau,step_interval):
		# print("Fitting %d, %d of %d" % (i, i+tau, len(x)));
		x_cut=x[i:i+tau];
		y_cut=y[i:i+tau];
		sig_cut=sig[i:i+tau];
		params, cov = linear_fitting_menke(x_cut, y_cut, np.mean(sig_cut), verbose=0);
		slopes.append(params[0]);		
	
	consec_differences=[];	
	for i in range(len(slopes)-1):
		consec_differences.append(slopes[i]-slopes[i+1]);
	AVR=0.5*np.var(consec_differences);
	uncertainty = np.sqrt(AVR);
	covm=[[AVR, 0], [0,0]];
	perr=[np.sqrt(AVR),0];
	if verbose:
		print("# Fit parameters and parameter errors from AVR method :")
		print("pfit = ", params_overall);
		print("perr = ", perr);

	return params_overall, covm;



