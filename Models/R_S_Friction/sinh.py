# What happens when you change things according to 
# v = 2*vo*sinh(dtau/((a-b)*sigma_eff))
# to solve for the new normal stress when steady state velocity drops? 
# p=deltatau/(a-b). we assume this doesn't change. 
# sigmaeff = effective normal stress on the interface (we assume this changes)

# if p is larger or sigma_eff is smaller, the velocity is more sensitive to stress changes. 

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc

def driver():

	[parray, sigma_eff, velred] = configure(); # get a range of velocities and friction conditions
	
	# Compute the arrays that we're going to plot

	stress_results=[];
	for i in range(len(parray)):
		tempstresses=[];
		for j in range(len(velred)):
			sigmastar=get_sigmastar(parray[i],sigma_eff,velred[j]);
			tempstresses.append(sigmastar);
		stress_results.append(tempstresses);

	make_plots(parray, velred, stress_results);
	return;

# Set up the problem
def configure():
	parray=[0.01,0.5, 1.0, 1.5,3.0]; # the different friction and stress conditions
	sigma_eff=1;
	velred=np.arange(0.3,1.8,0.1);
	return [parray, sigma_eff, velred]; 

# The compute functions
def get_sigmastar(p,sigma,velred):
	sigmastar=p/(isinh(velred*np.sinh(p/sigma) ));
	return sigmastar;

def isinh(x):
	isinh=np.log(x+np.sqrt(1+x*x));
	return isinh;



def make_plots(parray, velred, stress_results):
	plt.figure();
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	for i in range(len(parray)):
		plt.plot(stress_results[i],velred,label='p='+str(parray[i])+'MPa');
	plt.grid(True);
	plt.ylabel(r'Normalized Velocity',fontsize=16);
	plt.xlabel(r'Effective Normal Stress (MPa)',fontsize=16);
	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.legend();

	plt.text(1.54,1.04,r'$v=2v_0 \sinh{\frac{\Delta\tau}{(a-b)\bar{\sigma}}}$', fontsize=16, color='k');
	plt.savefig('velreductions.eps');
	return;

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})



if __name__=="__main__":
	driver();