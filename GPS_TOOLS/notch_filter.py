# Notch filter

# Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

import numpy as np
import matplotlib.pyplot as plt 
import collections
import copy


Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm


def notchfilt(x,fs,fn,Bn,filtfiltopt=True):
	# %
	# %NOTCHFILT kills one frequency content from the input signal
	# %   y = notchfilt(x,fs,fn,Bn,<filtfiltopt>) kills the frequency content fn 
	# %   of the signal x with the notch bandwidth of Bn.
	# %
	# %INPUT:
	# %   x       1-D signal array
	# %   fs      sampling frequency, Hz
	# %   fn      notch frequency, Hz
	# %   Bn      notch bandwidth, Hz
	# %   filtfiltopt   filtfilt option (to compensate phase distortion)
	# %                 default is true
	# %
	# %OUTPUT:
	# %   y       notch-filtered output signal
	# %   Sang-Ho Yun, 2007/07/27


	w0 = 2*np.pi*(fn/fs);
	Bw = 2*np.pi*(Bn/fs);

	alpha = (1-np.tan(Bw/2.0))/(1+np.tan(Bw/2.0));
	beta = np.cos(w0);

	a0 = 1;
	a1 = -beta*(1+alpha);
	a2 = alpha;

	b0 = (1+alpha)/2;
	b1 = -beta*(1+alpha);
	b2 = (1+alpha)/2;

	y = np.zeros(len(x));
	y[0] = x[0];
	y[1] = x[1];
	for n in range(2,len(x)):
		y[n] = 1.0/a0*(b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] -a2*y[n-2]);

	if filtfiltopt==True:
		y = y[::-1];  # reverse the direction. 
		x = copy.deepcopy(y);

		for n in range(2,len(x)):
			y[n] = 1.0/a0*(b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] -a2*y[n-2]);
		y = y[::-1];  # reverse the direction.

	return y;




def notch_filter_example():
	# The same as Paul/Sang-Ho's example for plotting the results of the notch filter. 
	# 
	N = 5e4;         # number of samples
	dt = 1;          # sampling interval (sec)
	t = np.arange(0,N,dt);
	fs = 1.0/dt;       # sampling frequency (Hz)
	df = fs/N;
	f = np.arange(-N/2.0, N/2.0, 1);
	f = [i*df for i in f];

	f0 = 1/(60*60);  # signal frequency (1cycle/hour)
	x = np.sin(2*np.pi*f0*t);
	X = np.abs(np.fft.fft(x,N));
	Xdc = X[0];

	Bw = 0.1*f0;

	y = notchfilt(x,fs,f0,Bw,True);
	Y = np.abs(np.fft.fft(y,N));
	Ydc = Y[0];


	plt.figure();
	plt.plot(t,x); plt.grid(True);
	plt.title('input signal'); 
	plt.xlabel('time (sec)'); 
	plt.ylabel('amplitude');
	plt.plot(t,y, 'r');
	plt.title('filtered signal'); 
	plt.xlabel('time (sec)'); 
	plt.ylabel('amplitude');
	plt.savefig('testnotch.eps');

	plt.figure(); 
	plt.plot(f,20*np.log10(np.fft.fftshift(X)/Xdc)); 
	plt.grid(True);
	plt.title('input and filtered spectra'); plt.xlabel('frequency (Hz)'); plt.ylabel('power (dB)');
	plt.xlim([0, 5*f0]);
	plt.plot(f,20*np.log10(np.fft.fftshift(Y)/Ydc), 'r'); 
	plt.savefig('notchspectra.eps');

	return;





