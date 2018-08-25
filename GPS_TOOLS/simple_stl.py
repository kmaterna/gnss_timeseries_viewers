# A set of functions to remove seasonal components from GPS time series using the STL algorithm
# Based on Cleveland et al., 1990
# And implementatoin by Daniel Trugman, 2014

import numpy as np 


def stl_simple(Y, n_p=365, n_s=41, n_i=2):

# Simplified Matlab Implementation for STL algorithm (Cleveland et al.
# 1990). Assumes uniformly spaced, complete (no gaps) time series.
# Daniel Trugman, 11/2014

#% Input Parameters
#    % Y: data vector - uniformly spaced column, complete (no gaps)
#    % n_p: number of epochs in a cycle
#        % optional, default is 365 (i.e., yearly cycle)
#    % n_s: seasonal smoothing span
#        % optional, default is 41
#    % n_i: inner loop iterations
#        % optional, default is 2

#% Output Parameters: Y = S + T + R
#    % S: seasonal component
#    % T: trend component
#    % R: residual component


#%%%%% derive STL parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# n_l:low-pass filter smoothing span 
if mod(n_p, 2) == 0:
	n_l = n_p;
else:
	n_l = n_p + 1;
# n_t: trend smoothing span
n_t = np.ceil(1.5*n_p/(1-1.5/n_s));
if mod(n_t, 2) == 0:
	n_t = n_t + 1;
# n_o: outer loop iterations
n_o = 0; # set to zero for now- no outer loop
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# format data
# Y = Y(:); # full time series as a column vector
n_n = len(Y); # total number of data entries
n_c = np.ceil(n_n/n_p); # number of cycles
    
# make sure have integer number of cycles
if mod(n_n/n_p,1) > 0.01:
	print('STL error! Incomplete final cycle');
	S = NaN(size(Y)); T = NaN(size(Y)); R = NaN(size(Y));
	return;
SC = np.reshape(Y, n_c, n_p);
sz_SC = size(SC);

# initialize seasonal, trend, and subcycle matrices
Sk = 0*Y; Tk = Sk;
SCk = 0*SC;
Lk = Sk; # low-pass filter vector

# inner loop: iteratively estimate S and T
for kk in range(0,n_i-1):  # matlab: 1:n_i
    
	# 1. Detrend the time series
	Sk = Y-Tk;
	SCk = np.reshape(Sk, sz_SC);
    
    # 2. Smooth the cycle subseries
	for jj in range(0,n_c-1): # matlab: 1:n_c
		SCk(jj,:) = smooth(SCk(jj,:), n_s, 'rloess');
	Sk = SCk(:);
    
    # 3. Low-pass filter the smoothed seasonal
	Lk = smooth(Sk, n_p); # moving avg. filter
	Lk = smooth(Lk, n_p); # another moving avg. filter
	Lk = smooth(Lk,3); # a third moving avg. filter with smaller window
	Lk = smooth(Lk, n_l, 'rloess'); # rloess smoothing with n_l window
    
    # 4. Detrend the smoothed seasonal
	Sk = Sk - Lk;
    
    # 5. Deseasonalize the trend
	Tk = Y - Sk;
    
    # 6. Smooth the trend
	Tk = smooth(Tk, n_t, 'rloess');
    
	print('\n');
	print('STL inner loop iteration ' +str(kk)+ ' complete.');

# finalize S, T, R
S = Sk; T = Tk;
R = Y - S - T;

print('STL decomposition complete.');

return [S, T, R];

