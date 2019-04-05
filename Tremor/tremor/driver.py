# tremor driver.py

import datetime as dt
import tremor_io
import tremor_tools
import tremor_plots


if __name__=="__main__":

	# A detailed time-tremor plot
	tremortype='wech_custom';
	tremor = tremor_tools.read_custom_tremor(tremortype);
	tremor_io.write_tremor_as_txt(tremor,"wech_2019_mod_catalog.txt");
