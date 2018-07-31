
import numpy as np
import collections
import datetime as dt 


Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);


def read_pbo_vel_file(infile):
# Meant for reading velocity files from the PBO/UNAVCO website. 
# Returns a Velfield collections object. 
	start=0;
	ifile=open(infile,'r');
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	for line in ifile:
		if start==1:
			temp=line.split();
			name.append(temp[0]);
			nlat.append(float(temp[7]));
			elon_temp=float(temp[8]);
			if elon_temp>180:
				elon_temp=elon_temp-360.0;
			elon.append(elon_temp);
			n.append(float(temp[19])*1000.0);
			e.append(float(temp[20])*1000.0);
			u.append(float(temp[21])*1000.0);
			sn.append(float(temp[22])*1000.0);
			se.append(float(temp[23])*1000.0);
			su.append(float(temp[24])*1000.0);
			t1=temp[-2];
			t2=temp[-1];
			first_epoch.append(int(t1[0:8]));
			last_epoch.append(int(t2[0:8]));
		if "*" in line:
			start=1;
	ifile.close();
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield];



def clean_velfield(velfield, num_years, max_sigma, coord_box):
# Take the raw GPS velocities, and clean them up. 
# Remove data that's less than num_years long, 
# has formal uncertainties above max_sigma, 
# or is outside our box of intersest. 
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	for i in range(len(velfield.n)):
		if velfield.sn[i] > max_sigma:  
			continue;
		if velfield.se[i] > max_sigma:
			continue;
		if velfield.last_epoch[i]-velfield.first_epoch[i] <= num_years*10000:
			continue;
		if (velfield.elon[i]>coord_box[0] and velfield.elon[i]<coord_box[1] and velfield.nlat[i]>coord_box[2] and velfield.nlat[i]<coord_box[3]):
			#The station is within the box of interest. 
			name.append(velfield.name[i]);
			nlat.append(velfield.nlat[i]);
			elon.append(velfield.elon[i]);
			n.append(velfield.n[i]);
			sn.append(velfield.sn[i]);
			e.append(velfield.e[i]);
			se.append(velfield.se[i]);
			u.append(velfield.u[i]);
			su.append(velfield.su[i]);			
			first_epoch.append(velfield.first_epoch[i]);
			last_epoch.append(velfield.last_epoch[i]);
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield];


def remove_duplicates(velfield):
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	
	for i in range(len(velfield.n)):
		is_duplicate = 0;
		for j in range(len(name)):
			if abs(nlat[j]-velfield.nlat[i])<0.0005 and abs(elon[j]-velfield.elon[i])<0.0005:
				# we found a duplicate measurement. 
				is_duplicate = 1;
				# Right now assuming all entries at the same lat/lon have the same velocity values. 

		if is_duplicate == 0:
			name.append(velfield.name[i]);
			nlat.append(velfield.nlat[i]);
			elon.append(velfield.elon[i]);
			n.append(velfield.n[i]);
			sn.append(velfield.sn[i]);
			e.append(velfield.e[i]);
			se.append(velfield.se[i]);
			u.append(velfield.u[i]);
			su.append(velfield.su[i]);			
			first_epoch.append(velfield.first_epoch[i]);
			last_epoch.append(velfield.last_epoch[i]);

	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);	
	return [myVelfield];



def read_pbo_pos_file(filename):
	[yyyymmdd, Nlat, Elong, dN, dE, dU, Sn, Se, Su] = np.loadtxt(filename, skiprows=37, unpack=True,usecols=(0,12,13,15,16,17,18,19,20));
	dN=[i*1000.0 for i in dN];
	dE=[i*1000.0 for i in dE];
	dU=[i*1000.0 for i in dU];
	Sn=[i*1000.0 for i in Sn];
	Se=[i*1000.0 for i in Se];
	Su=[i*1000.0 for i in Su];
	specific_file=filename.split('/')[-1];
	dtarray= [dt.datetime.strptime(str(int(i)),"%Y%m%d") for i in yyyymmdd];
	myData=Timeseries(name=specific_file[0:4],coords=[Elong[0]-360, Nlat[0]], dtarray=dtarray, dN=dN, dE=dE, dU=dU, Sn=Sn, Se=Se, Su=Su, EQtimes=[]);
	return [myData];



def read_noel_file(filename):
	
	names = np.genfromtxt(filename,skip_header=8, usecols=(0), dtype=('unicode') );
	E, N, U, Ea1, Na1, Ua1, Ea2, Na2, Ua2, Es1, Ns1, Us1, Es2, Ns2, Us2 = np.genfromtxt(filename,skip_header=8, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True );

	return [names, E, N, U, Ea1, Na1, Ua1, Ea2, Na2, Ua2, Es1, Ns1, Us1, Es2, Ns2, Us2];


def read_noel_file_station(filename,station):
	names = np.genfromtxt(filename,skip_header=8, usecols=(0), dtype=('unicode') );
	E, N, U, Ea1, Na1, Ua1, Ea2, Na2, Ua2, Es1, Ns1, Us1, Es2, Ns2, Us2 = np.genfromtxt(filename,skip_header=8, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True );
	i=np.where(names==station)[0][0];
	return [E[i], N[i], U[i], Ea1[i], Na1[i], Ua1[i], Ea2[i], Na2[i], Ua2[i], Es1[i], Ns1[i], Us1[i], Es2[i], Ns2[i], Us2[i]];

