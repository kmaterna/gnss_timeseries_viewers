# tremor driver.py

import datetime as dt
import tremor_io
import tremor_tools
import tremor_plots


if __name__=="__main__":
	readfuncs={"wech":tremor_io.read_wech,
	"ide":tremor_io.read_ide};
	filenames={"wech":"../../GPS_POS_DATA/tremor/08_01_2009_10_31_2018.txt",
	"ide":"../../GPS_POS_DATA/tremor/trm_Cascadia.20050101.3652.92921871.csv"};
	#"wech":"../../GPS_POS_DATA/tremor/08_01_2009_10_31_2018.txt",

	# A detailed time-tremor plot
	tremortype='wech';
	tremor=readfuncs[tremortype](filenames[tremortype]);
	tremor_plots.timing_2014_plot(tremor, tremortype);


	# Print the tremor by time, in hours or days from event time, used for GMT
	# Really this should not be in this section of code, this is very disorganized
	ETS14time=dt.datetime.strptime('20140210 00',"%Y%m%d %H");  # time of 2014 earthquake
	eqtime14=dt.datetime.strptime('20140310 05',"%Y%m%d %H");  # time of 2014 earthquake
	eqtime16=dt.datetime.strptime('20161208 00',"%Y%m%d %H");  # time of 2016 earthquake
	tremor_tools.write_deltat_certain_day(tremor,ETS14time,ETS14time+dt.timedelta(days=27),"Feb2014.txt",'days');
	tremor_tools.write_deltat_certain_day(tremor,eqtime14,eqtime14+dt.timedelta(days=5),"March10.txt",'hours');
	tremor_tools.write_deltat_certain_day(tremor,eqtime16,eqtime16+dt.timedelta(days=100),"Dec08.txt",'days');
