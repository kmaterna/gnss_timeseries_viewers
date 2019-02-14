# tremor driver.py

import datetime as dt
import tremor_io
import tremor_tools
import tremor_plots


if __name__=="__main__":

	# A detailed time-tremor plot
	tremortype='wech';
	tremor_io.read_input_tremor(tremortype);
	tremor_plots.timing_2014_plot(tremor, tremortype);

	# Print the tremor by time, in hours or days from event time, used for GMT
	# Really this should not be in this section of code, this is very disorganized
	ETS14time=dt.datetime.strptime('20140210 00',"%Y%m%d %H");  # time of 2014 earthquake
	eqtime14=dt.datetime.strptime('20140310 05',"%Y%m%d %H");  # time of 2014 earthquake
	eqtime16=dt.datetime.strptime('20161208 00',"%Y%m%d %H");  # time of 2016 earthquake
	tremor_tools.write_deltat_certain_day(tremor,ETS14time,ETS14time+dt.timedelta(days=27),"Feb2014.txt",'days');
	tremor_tools.write_deltat_certain_day(tremor,eqtime14,eqtime14+dt.timedelta(days=5),"March10.txt",'hours');
	tremor_tools.write_deltat_certain_day(tremor,eqtime16,eqtime16+dt.timedelta(days=100),"Dec08.txt",'days');
