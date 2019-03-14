import numpy as np 
import matplotlib.pyplot as plt
import collections
import datetime as dt 
import tremor_io
import tremor_tools

TremorCat = collections.namedtuple("TremorCat",['dtarray','lonarray','latarray']);

def configure():
	tremor_type = "wech_custom_no_depth";
	bounds = [-125,-120,38,43];
	spacing = 0.1; 
	outfile="tremor_density.txt";
	return [tremor_type, bounds, spacing, outfile];

def compute(tremor, bounds, spacing):
	xarray = np.arange(bounds[0], bounds[1], spacing);
	yarray = np.arange(bounds[2], bounds[3], spacing);
	print(yarray);
	print(xarray);
	density = np.zeros([len(yarray), len(xarray)]);
	for i in range(len(yarray)):
		for j in range(len(xarray)):
			box_interest = [xarray[j], xarray[j]+spacing, yarray[i], yarray[i]+spacing];
			restricted = tremor_tools.restrict_to_box(tremor, box_interest);
			density[i][j] = len(restricted.lonarray);
	return xarray, yarray, density;

def outputs(xarray, yarray, spacing, density, outfile):
	ofile=open(outfile,'w');
	for i in range(len(yarray)):
		for j in range(len(xarray)):
			xpos = xarray[j]+spacing/2.0;
			ypos = yarray[i]+spacing/2.0;
			ofile.write("%f %f %d\n" % (xpos, ypos, density[i][j]) );
	ofile.close();
	return;

def plot_outfile(outfile):
	[lon, lat, density] = np.loadtxt(outfile,unpack=True);
	plt.figure(figsize=(8,6));
	plt.scatter(lon, lat, c=density, s=60, marker='s');
	plt.colorbar();
	plt.savefig('tremor_density');
	return;

if __name__=="__main__":
	[tremor_type, bounds, spacing, outfile] = configure();
	# tremor = tremor_tools.read_custom_tremor(tremor_type);
	# tremor = TremorCat(dtarray=[dt.datetime.strptime("20140202","%Y%m%d")], lonarray=[-122.9],latarray=[40.3]); # a test case 
	# xarray, yarray, density = compute(tremor, bounds, spacing);
	# outputs(xarray, yarray, spacing, density, outfile);
	plot_outfile(outfile);
