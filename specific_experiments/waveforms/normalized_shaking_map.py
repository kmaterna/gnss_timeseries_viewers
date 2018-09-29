import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.imaging.beachball import beach
import sys
from mpl_toolkits.basemap import Basemap
import haversine
import shaking_map


def run_map(norm_eqcode, exp_eqcode):
	[lonbounds, latbounds, normfile,expfile,new_ampfile,eqfile]=configure(norm_eqcode, exp_eqcode);
	[norm_namecodes, norm_lons, norm_lats, norm_radials, norm_transverses, norm_verticals]=shaking_map.inputs_amplitudes(normfile);
	[exp_namecodes, exp_lons, exp_lats, exp_radials, exp_transverses, exp_verticals]=shaking_map.inputs_amplitudes(expfile);
	[new_namecodes, new_lons, new_lats, new_radials, new_transverses, new_verticals]=compute(norm_namecodes, norm_lons, norm_lats, norm_radials, norm_transverses, norm_verticals,exp_namecodes, exp_lons, exp_lats, exp_radials, exp_transverses, exp_verticals);
	shaking_map.write_amplitudes(new_namecodes, new_lons, new_lats, new_radials, new_transverses, new_verticals, new_ampfile);
	shaking_map.amplitude_map(eqfile, new_ampfile, lonbounds, latbounds,exp_eqcode,'verticals','Amplitude_maps/norm_'+exp_eqcode); # THE BASEMAP WITH STATIONS
	shaking_map.amplitude_map(eqfile, new_ampfile, lonbounds, latbounds,exp_eqcode,'radials','Amplitude_maps/norm_'+exp_eqcode); # THE BASEMAP WITH STATIONS
	shaking_map.amplitude_map(eqfile, new_ampfile, lonbounds, latbounds,exp_eqcode,'tangentials','Amplitude_maps/norm_'+exp_eqcode); # THE BASEMAP WITH STATIONS
	return;

def configure(norm_eqcode, exp_eqcode):
	latbounds=[38.8,42];
	lonbounds=[-126.5,-121];
	normfile="Amplitude_maps/ALL/outputs_"+norm_eqcode+".txt";
	expfile="Amplitude_maps/ALL/outputs_"+exp_eqcode+".txt";
	new_ampfile="Amplitude_maps/ALL/normalized_"+exp_eqcode+".txt"
	eqfile="Amplitude_maps/eq_database.txt";
	return [lonbounds, latbounds, normfile,expfile,new_ampfile,eqfile];

def compute(norm_namecodes, norm_lons, norm_lats, norm_radials, norm_transverses, norm_verticals,exp_namecodes, exp_lons, exp_lats, exp_radials, exp_transverses, exp_verticals):
	new_namecodes=[]; new_lons=[]; new_lats=[]; new_radials=[]; new_transverses=[]; new_verticals=[];

	for i in range(len(exp_namecodes)):
		if exp_namecodes[i] in norm_namecodes:
			station_index=norm_namecodes.index(exp_namecodes[i]);
			normalizing_value=norm_verticals[station_index];
			if ~np.isnan(normalizing_value) and normalizing_value>0.0:
				new_namecodes.append(exp_namecodes[i]);
				new_lons.append(exp_lons[i]);
				new_lats.append(exp_lats[i]);
				new_radials.append(exp_radials[i]/normalizing_value);
				new_transverses.append(exp_transverses[i]/normalizing_value);
				new_verticals.append(exp_verticals[i]/normalizing_value);
	return [new_namecodes, new_lons, new_lats, new_radials, new_transverses, new_verticals];

if __name__=="__main__":
	# run_map("2017","2016");
	run_map("2017","2014");
	# run_map("2017","2010");
	# run_map("2017","2005");

