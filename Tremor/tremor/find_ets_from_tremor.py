# Identify jumps in the cumulative tremor distribution


import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import tremor_io
import tremor_tools

def find_ets(tremor):
	start_time=dt.datetime.strptime('20120301',"%Y%m%d");
	end_time=dt.datetime.strptime('20181101',"%Y%m%d");
	box_interest=[-123.3,-123,40,41];

	tremor=tremor_tools.restrict_to_box(tremor, box_interest, start_time, end_time);
	[dts, rate] = tremor_tools.get_rates(tremor);

	for i in range(len(dts)):
		if rate[i]>20:
			print(dts[i]);

	plt.figure();
	plt.plot_date(dts,rate,linestyle='-',marker=None,color='b',linewidth=4);
	plt.xlabel('Time',fontsize=20);
	plt.ylabel('Daily Tremor Rates',fontsize=20);
	plt.savefig('rates.eps');

	return;






if __name__=="__main__":
	tremor=tremor_io.read_wech("../../GPS_POS_DATA/tremor/08_01_2009_10_31_2018.txt");
	find_ets(tremor);