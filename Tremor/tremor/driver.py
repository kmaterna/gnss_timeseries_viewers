# tremor driver.py

import datetime as dt
import tremor_io
import tremor_tools
import tremor_plots


if __name__=="__main__":

	# A detailed time-tremor plot
	tremortype='wech_custom';
	tremor = tremor_io.read_input_tremor(tremortype);
	tremor_plots.simple_plot_2015_2017(tremor, tremortype);

	# tremor = tremor_tools.read_custom_tremor(tremortype);  # this is the combined wech + wech_custom catalog
	# tremor_io.write_tremor_as_txt(tremor,"wech_2019_mod_catalog.txt");
