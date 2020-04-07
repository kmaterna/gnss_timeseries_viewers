# Project Rainbow Vector
# Plot a vector that changes as a function of time
# Along with the time series that it came from. 
# For viewing the evolution of a transient 
# April 2020


import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import datetime as dt 




def plotting_function(dataobj_list, sorted_distances, params):
	# Setting annotations
	# # This is good for SoCal
	EQtimes=[];
	EQtimes.append(dt.datetime.strptime("20100403", "%Y%m%d"));  # starts with the most important one
	EQtimes.append(dt.datetime.strptime("20120826", "%Y%m%d"));	
	start_time_plot=dt.datetime.strptime("20050101","%Y%m%d");
	end_time_plot=dt.datetime.strptime("20200116", "%Y%m%d");
	label_date=dt.datetime.strptime("20200215","%Y%m%d");
	offset=0;
	spacing=20;	
	vmax=5000;
	chosen_colormap='jet_r' # jet_r is okay. Also tried gnuplot, rainbow, and brg
	color_boundary_object=matplotlib.colors.Normalize(vmin=0,vmax=vmax, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap=chosen_colormap);  # maybe useful? 

	fig, axarr = plt.subplots(2,2,figsize=(15,13),dpi=300);

	# East Plots
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		edata=dataobj_list[i].dE;
		emean=np.nanmean(dataobj_list[i].dE);
		edata=[x + offset - emean for x in edata];
		date_deltas=[];
		for j in range(len(edata)):
			deltas = dataobj_list[i].dtarray[j] - start_time_plot;
			date_deltas.append(deltas.days);
		l1 = axarr[0][0].scatter(dataobj_list[i].dtarray,edata,c=date_deltas,marker='+',s=1.5,cmap=chosen_colormap,vmin=0, vmax=vmax);
		axarr[0][0].text(label_date,edata[0],dataobj_list[i].name,fontsize=9,color='black');
	axarr[0][0].set_xlim(start_time_plot,end_time_plot);
	axarr[0][0].set_ylim([-20,offset+15])
	bottom,top=axarr[0][0].get_ylim();
	for i in range(len(EQtimes)):
		axarr[0][0].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	axarr[0][0].set_ylabel("East (mm)");
	axarr[0][0].set_title("East GPS Time Series")
	axarr[0][0].grid(True)


	# North Plots
	offset=0;
	for i in range(len(dataobj_list)):
		offset=spacing*i;
		ndata=dataobj_list[i].dN;
		nmean=np.nanmean(dataobj_list[i].dN);
		ndata=[x + offset - nmean for x in ndata];
		date_deltas=[];
		for j in range(len(ndata)):
			deltas = dataobj_list[i].dtarray[j] - start_time_plot;
			date_deltas.append(deltas.days);
		l1 = axarr[0][1].scatter(dataobj_list[i].dtarray,ndata,c=date_deltas,marker='+',s=1.5,cmap=chosen_colormap, vmin=0, vmax=5000);
	axarr[0][1].set_xlim(start_time_plot,end_time_plot);
	axarr[0][1].set_ylim([-20,offset+15])
	bottom,top=axarr[0][1].get_ylim();
	for i in range(len(EQtimes)):
		axarr[0][1].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	axarr[0][1].set_ylabel("North (mm)");
	axarr[0][1].set_title("North GPS Time Series")
	axarr[0][1].grid(True)

	# Vert Plots
	offset=0;
	for i in range(len(dataobj_list)):
		offset=spacing*2*i;
		udata=dataobj_list[i].dU;
		umean=np.nanmean(dataobj_list[i].dU);
		udata=[x + offset - umean for x in udata];
		date_deltas=[];
		for j in range(len(udata)):
			deltas = dataobj_list[i].dtarray[j] - start_time_plot;
			date_deltas.append(deltas.days);
		l1 = axarr[1][0].scatter(dataobj_list[i].dtarray,udata,c=date_deltas,marker='+',s=1.5,cmap=chosen_colormap, vmin=0, vmax=5000);
		axarr[1][0].text(label_date,udata[0],dataobj_list[i].name,fontsize=9,color='black');
	axarr[1][0].set_xlim(start_time_plot,end_time_plot);
	axarr[1][0].set_ylim([-20,offset+15])
	bottom,top=axarr[1][0].get_ylim();
	for i in range(len(EQtimes)):
		axarr[1][0].plot_date([EQtimes[i], EQtimes[i]], [bottom, top], '--k'); 
	axarr[1][0].set_ylabel("Vertical (mm)");
	axarr[1][0].set_title("Vertical GPS Time Series")
	axarr[1][0].grid(True)

	# Vector Map
	for i in range(len(dataobj_list)):
		axarr[1][1].plot(dataobj_list[i].coords[0], dataobj_list[i].coords[1],marker='v',markersize=6,color='black');
		[evec, nvec, uvec, num_days] = make_vector(dataobj_list[i],start_time_plot);
		for j in range(len(evec)):
			mycolor = custom_cmap.to_rgba(num_days[j]);
			axarr[1][1].quiver(dataobj_list[i].coords[0], dataobj_list[i].coords[1], evec[j], nvec[j], color=mycolor,scale=60);

	axarr[1][1].quiver(-115.75, 32.79, 10, 0, color=custom_cmap.to_rgba(3600),scale=60);
	axarr[1][1].text(-115.75, 32.80, "10mm");
	axarr[1][1].plot(-115.5153,33.0153,marker='*',markersize=14, color='red');
	axarr[1][1].text(-115.50,33.0153,"NBGF",fontsize=14);
	axarr[1][1].set_ylim([32.75, 33.35]);
	axarr[1][1].set_xlim([-115.8, -115.2]);

	plt.savefig(params.outdir+"/"+params.outname+'_RainbowVector.png');
	return;



def make_vector(Data0, start_time):
	# Every once in a while, extract a date and a vector displacement
	# Return those values. 
	# Very simple right now, no smoothing or filtering. 
	evec = []; 
	nvec = [];
	uvec = [];
	num_days = [];
	sampling_interval = 200;
	inc=45;

	for i in range(len(Data0.dtarray)):
		if Data0.dtarray[i]>dt.datetime.strptime("20180101","%Y%m%d"):
			break;
		date_delta = Data0.dtarray[i]-start_time;
		if date_delta.days % sampling_interval == 0 :

			num_days.append(date_delta.days);
			evec.append(np.nanmedian(Data0.dE[i-inc:i+inc]));
			nvec.append(np.nanmedian(Data0.dN[i-inc:i+inc]));
			uvec.append(np.nanmedian(Data0.dU[i-inc:i+inc]));

	return [evec, nvec, uvec, num_days];




