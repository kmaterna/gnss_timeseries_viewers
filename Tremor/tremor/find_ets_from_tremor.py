# Identify jumps in the cumulative tremor distribution

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import tremor_io
import tremor_tools

def find_ets(tremor):
	start_time=dt.datetime.strptime('20120301',"%Y%m%d");
	end_time=dt.datetime.strptime('20181101',"%Y%m%d");
	box_interest=[-123.3,-123,40.2,40.8];

	tremor=tremor_tools.restrict_to_box(tremor, box_interest, start_time, end_time);
	[dts, rate] = tremor_tools.get_rates(tremor);

	ets_dates=[];
	for i in range(len(dts)):
		if rate[i]>9:
			print(dts[i]);
			ets_dates.append(dts[i]);
	
	# Decluster	
	print("Declustering the ETS days");
	ets_declustered=[];	
	post_period=[];

	for i in range(len(ets_dates)-1):
		post_dt=ets_dates[i+1]-ets_dates[i];
		post_days=post_dt.days;
		if post_days>60:  # ignoring the little noise
			ets_declustered.append(ets_dates[i]);
			post_period.append(post_days);
			print(ets_dates[i]);
			print("Followed by %d days of quiet" % post_days);

	print("Average post-period plus or minus:");
	print(np.mean(post_period));
	print(np.std(post_period));

	plt.figure();
	plt.plot_date(dts,rate,linestyle='-',marker=None,color='b',linewidth=4);
	plt.xlabel('Time',fontsize=20);
	plt.ylabel('Daily Tremor Rates',fontsize=20);
	plt.grid(True);
	plt.savefig('rates.eps');

	return;




if __name__=="__main__":
	tremor = tremor_tools.read_custom_tremor("wech_custom");
	find_ets(tremor);
