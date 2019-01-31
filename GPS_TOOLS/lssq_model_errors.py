# A toolbox for fitting by least squares and getting uncertainties on velocity parameters. 
# TOOLS:
# Linear_fitting_menke  : the formal uncertainties from the textbook
# fit_curvefit          : a python function that does practically the same thing, but with more sophisticated handling of individual errors
# colored_noise         : something I'm writing now. 

import numpy as np 
from scipy import optimize


def linear_fitting_menke(x, y, sig):
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
	print("pfit = ", params);
	print("perr = ", perr)

	return params, covm;





def fit_curvefit(x, y, sig):
	"""
	Note: As per the current documentation (Scipy V1.1.0), sigma (yerr) must be:
	None or M-length sequence or MxM array, optional
	p0 is an initial guess for the first parameter
	f is the function we're fitting (line in this case)
	sigma must be an array
	"""

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
	
	print("# Fit parameters and parameter errors from curve_fit method :")
	print("pfit = ", pfit);
	print("perr = ", perr);

	return pfit, pcov;



 
	