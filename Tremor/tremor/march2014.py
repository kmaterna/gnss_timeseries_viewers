# Make maps of tremor on 3/10/2014

import numpy as np
import datetime as dt  
import tremor_io


def see_march_2014(tremor):
	ofile=open('march.txt','w');
	for i in range(len(tremor.dtarray)):
		if tremor.dtarray[i]>=dt.datetime.strptime("20140310","%Y%m%d") and tremor.dtarray[i]<=dt.datetime.strptime("20140313","%Y%m%d"):
			ofile.write("%s %f %f\n" % (dt.datetime.strftime(tremor.dtarray[i],"%Y%m%d"),tremor.lonarray[i], tremor.latarray[i]) );
			print(tremor.lonarray[i])
	ofile.close();
	return;



if __name__=="__main__":
	readfuncs={"wech":tremor_io.read_wech,
	"ide":tremor_io.read_ide};
	filenames={"wech":"../../GPS_POS_DATA/tremor/08_01_2009_10_31_2018.txt",
	"ide":"../../GPS_POS_DATA/tremor/trm_Cascadia.20050101.3652.92921871.csv"};

	tremortype='wech'
	tremor=readfuncs[tremortype](filenames[tremortype]);

	see_march_2014(tremor);