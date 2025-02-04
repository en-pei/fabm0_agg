import matplotlib
from matplotlib import rc
from matplotlib import ticker, cm
import cftime
import netCDF4
from pylab import *
from scipy.interpolate import interp1d
#import netcdftime
import numpy	
import os
import sys
import subprocess
from statistics import mean
import re
import pickle
#from glob import glob

####code for postprocessing and plotting of multiple runs in the same parent folder#####

filelist=[]
path='/gpfs/work/lie/run/par/4param'  #/run/par/test
#path='/gpfs/work/lie/run/par/6oszmay' #6param'  #linked best run from parameterauto.sh

#path='/gpfs/work/lie/run/par/oszmayruns'  #partest'
#path='/gpfs/home/lie/setups/bale2002/par'
#input from runpost.py to iterate through all layers
lay=sys.argv[1] #5 #16 #4 #16 #6 #16 #layer for SPM concentration
laysize=sys.argv[2]# 5 #16 #14 #9 #layer for ESD
layflow=sys.argv[3]
#print(lay,laysize,'lay')

h=3600./60 #/120. output was per 120s
tstart=int(0*h)#int(6.*h)
tend=tstart+int(39*h)#tstart+int(8*h)
tslice=slice(tstart,tend)

datadir = '/gpfs/home/lie/setups/bale2002/Plot/data'
#data = asarray(numpy.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=', '))
#sizedata = asarray(numpy.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))
#s1datatime=sizedata[:,0]-5#-0.25
#s1datasize=sizedata[:,1] #data is in micrometer
#datatime=data[:,0]-5.#-0.25
#dataspmc=data[:,2]

#winteresd=asarray(numpy.loadtxt('%s/esd_feb_data.txt'%datadir,delimiter=', '))
#summeresd=asarray(numpy.loadtxt('%s/esd_july_data.txt'%datadir,delimiter=', '))
#winterspmc=asarray(numpy.loadtxt('%s/spmc_feb_data.txt'%datadir,delimiter=', '))
#summerspmc=asarray(numpy.loadtxt('%s/spmc_july_data.txt'%datadir,delimiter=', '))

#winteresdm=asarray(numpy.loadtxt('%s/esd_feb_model.txt'%datadir,delimiter=', '))
#summeresdm=asarray(numpy.loadtxt('%s/esd_july_model.txt'%datadir,delimiter=', '))
#winterspmcm=asarray(numpy.loadtxt('%s/spmc_feb_model.txt'%datadir,delimiter=', '))
#summerspmcm=asarray(numpy.loadtxt('%s/spmc_july_model.txt'%datadir,delimiter=', '))

bloomspmc=asarray(numpy.loadtxt('%s/bloomspm.txt'%datadir,delimiter=','))
regularspmc=asarray(numpy.loadtxt('%s/regularspm.txt'%datadir,delimiter=','))
bloomd50=asarray(numpy.loadtxt('%s/bloomd50.txt'%datadir,delimiter=','))
regulard50=asarray(numpy.loadtxt('%s/regulard50.txt'%datadir,delimiter=','))

tmzspmjuly=asarray(numpy.loadtxt('%s/tmzspmjuly.txt'%datadir,delimiter=','))
tmzd50july=asarray(numpy.loadtxt('%s/tmzd50july.txt'%datadir,delimiter=','))
tmzflowjuly=asarray(numpy.loadtxt('%s/tmzflowjuly.txt'%datadir,delimiter=','))
tmzspmfeb=asarray(numpy.loadtxt('%s/tmzspmfeb.txt'%datadir,delimiter=','))
tmzd50feb=asarray(numpy.loadtxt('%s/tmzd50feb.txt'%datadir,delimiter=','))
tmzflowfeb=asarray(numpy.loadtxt('%s/tmzflowfeb.txt'%datadir,delimiter=','))

oszspmmay=asarray(numpy.loadtxt('%s/oszspmmay.txt'%datadir,delimiter=','))
oszd50may=asarray(numpy.loadtxt('%s/oszd50may.txt'%datadir,delimiter=','))
oszflowmay=asarray(numpy.loadtxt('%s/oszflowmay.txt'%datadir,delimiter=','))
oszspmfeb=asarray(numpy.loadtxt('%s/oszspmfeb.txt'%datadir,delimiter=','))
oszd50feb=asarray(numpy.loadtxt('%s/oszd50feb.txt'%datadir,delimiter=','))
oszflowfeb=asarray(numpy.loadtxt('%s/oszflowfeb.txt'%datadir,delimiter=','))

bale3500spmc=asarray(numpy.loadtxt('%s/spmc_bale3500_new.txt'%datadir,delimiter=','))  #', '
bale3500spmcf=asarray(numpy.loadtxt('%s/flow_balespmc3500_new.txt'%datadir,delimiter=','))
bale3500d50=asarray(numpy.loadtxt('%s/d50_bale3500_new.txt'%datadir,delimiter=','))
bale3500d50f=asarray(numpy.loadtxt('%s/flow_baled503500_new.txt'%datadir,delimiter=','))#this is old


##define dataset
tsmdatat=bale3500spmc[8:,0]/20*16+4 -0.95-0.1-0.6 # oszspmmay[:,0]+6-0.3 #tmzspmjuly[:,0]+7.5 # bale3500spmc[:,0]-0.95-0.1-0.6 # +7 #tmzspmjuly[:,0]+7.5 #-0.42 #datatime #oszspmfeb[:,0] #tmzspmfeb[:-1,0]+4+6 #

d50datat=bale3500d50[:,0]/20*16+4 -0.95-0.1-0.6 #oszd50may[:,0]+6-0.3 #tmzd50july[:,0]+7.5 #bale3500d50[8:,0]-0.95-0.1-0.6 #+7 #tmzd50july[:,0]+7.5 #-0.42 #s1datatime #oszd50feb[:,0] #tmzd50feb[:,0]+4+6 #

flowdatat=bale3500spmcf[:,0]/20*16+4 -0.95-0.1-0.6 #oszflowmay[:,0]+6-0.3 #tmzflowjuly[:,0]+7.5 #bale3500d50f[8:,0]-0.95-0.1-0.6# +7 #tmzflowjuly[:,0]+7.5 #-0.42 #oszflowfeb[:,0] #tmzflowfeb[:,0]+4+6 #
#flowdatat1=bale3500spmcf[:,0]-0.95-0.1 #-0.42

tsmdata=bale3500spmc[0:,1] #oszspmmay[:,1] #tmzspmjuly[:,1] ##bale3500spmc[8:,1] #tmzspmjuly[:,1] #dataspmc #oszspmfeb[:,1] #tmzspmfeb[:-1,1] #
d50data=bale3500d50[0:,1] #oszd50may[:,1] #tmzd50july[:,1] ##bale3500d50[8:,1] #tmzd50july[:,1] #s1datasize #oszd50feb[:,1]  #tmzd50feb[:,1] #
flowdata=bale3500spmcf[0:,1] #oszflowmay[:,1] #tmzflowjuly[:,1] ##bale3500d50f[8:,1] #tmzflowjuly[:,1]#oszflowfeb[:,1] #tmzflowfeb[:,1] #

#tsmdatamax=bale3500spmc[0:,2]
#tsmdatamin=bale3500spmc[0:,3]
#d50datamax=bale3500d50[0:,2]
#d50datamin=bale3500d50[0:,3]
#flowdatamax=bale3500spmcf[0:,2]
#flowdatamin=bale3500spmcf[0:,3]

#flowdata1=bale3500spmcf[:,1] 

#s1datatimecorr=list(i for i in s1datatime) #(i-0.27) #time in datasize, 2.7 h earlier
#print len(datatime) #79
length=range(len(tsmdatat)) #datatime	len(winterspmc[:,0])			#idexing spmc data
lengths=range(len(d50datat)) 	#s1datasize len(winteresd[:,0])			#indexing size data
lengthf=range(len(flowdatat)) 	#s1datasize len(winteresd[:,0])			#indexing size data


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx#, array[idx]

def find_max(array):
    #Array = array#np.asarray(array)
    #maxidx = np.where(Array==max(array))
    maxidx=np.argmax(array)
    return maxidx 

Gidx=[]
Gsidx=[]
Gsidf=[]
filenumber=0
matplotlib.rcParams.update({'font.size': 22})
fig, axs = plt.subplots(6, sharex=True)  #4 (5, sharex=False)
rcParams['lines.linewidth']=5 #3
rcParams['lines.color']='firebrick' #3
for i in os.scandir(path):
	i.is_dir()
	filelist.append(i.name)
	try:
		ncfile=i.path+'/bale2002.nc'
		nc=netCDF4.Dataset(ncfile)
		ncv=nc.variables
		#utime=netcdftime.utime(ncv['time'].units)
		utime=cftime.utime(ncv['time'].units)
		time=utime.num2date(ncv['time'][:])
		secs=ncv['time'][:]
		#print len(secs)
		if len(secs)<481:  #here! for 4h runs!
			continue #skip run time smaller ones
		lpm=squeeze(ncv['spm_spm'][:,lay])
		agglpm = squeeze(ncv['agg_agglpm'][:,lay]) #[:,laysize]
#		phy=squeeze(ncv['npzd_phy'][:,lay])
#		det=squeeze(ncv['npzd_det'][:,lay])
		totlpm=lpm+agglpm#+phy+det+aggorg #no agglpm intended
		lpmsize = squeeze(ncv['spm_spm'][:,laysize])  #in size layer
		agglpmsize = squeeze(ncv['agg_agglpm'][:,laysize])
		totlpmsize=lpmsize+agglpmsize
		ctot=totlpm #squeeze(ncv['agg_Ctot'][:,lay]) 
		ctot1=totlpm #squeeze(ncv['agg_Ctot'][:,16]) 
		dsize=squeeze(ncv['agg_Dsize'][:,laysize]) #  dynamical size
		esd=dsize #squeeze(ncv['agg_esd'][:,laysize])  #dsize #dsize is for dynamical
#		esd1=squeeze(ncv['agg_Dresus'][:,laysize])#corrected #squeeze(ncv['agg_Dsize'][:,16]) #dynamical size
		Gsize = squeeze(ncv['agg_G'][:,laysize])
		G = squeeze(ncv['agg_G'][:,lay])
		ws = squeeze(ncv['agg_ws'][:,laysize]) #'agg_sinkinglpm' 'agg_resuspensionlpm'agg_diffsetlpm  (ncv['agg_ws'][:,lay])
		coagulation=squeeze(ncv['agg_coagulationlpm'][:,laysize])
		breakup=squeeze(ncv['agg_Breakup'][:,laysize])
		sinking=squeeze(ncv['agg_sinkinglpm'][:,laysize])
		diffset=squeeze(ncv['agg_diffsetlpm'][:,laysize])
		dD=squeeze(ncv['agg_resuspensionlpm'][:,laysize])
		xsize=squeeze(ncv['agg_Xsize'][:,laysize])

#		Dbreakup=squeeze(ncv['agg_Dbreakup'][:,laysize])
#		Dcoagulation=squeeze(ncv['agg_Dcoagulation'][:,laysize])
#		Dsinking=squeeze(ncv['agg_sinkinglpm'][:,laysize])
#		Derosion=squeeze(ncv['agg_dflux'][:])

		uu = abs(squeeze(ncv['u'][:,layflow]))*sqrt(2) #flow velocity

		#ws2 = squeeze(ncv['agg_sinkinglpm'][:,lay]) #agg_Breakup
		#ws = -ws1 +ws2
		aindex = np.diff(np.sign(np.diff(esd))).nonzero()[0] + 1             # local min & max
		ymin = (np.diff(np.sign(np.diff(esd))) > 0).nonzero()[0] + 1         # local min
		ymax = (np.diff(np.sign(np.diff(esd))) < 0).nonzero()[0] + 1         # local max# +1 
#		axs[0].plot(secs/3600,ctot,label=i.name) #,color='seagreen') #color='firebrick', totlpm    #+0.27 #SPMC	vs t
		axs[0].plot(secs/3600,coagulation,label=i.name)
		#axs[0].plot(secs/3600,phy+det,label=i.name) #totlpm    #+0.27 #SPMC	vs t
		#axs[0].plot(secs/3600,ctot1) #totlpm    #+0.27 #SPMC	vs t

#		axs[1].plot(secs/3600,esd*1e6,label=i.name) #,color='seagreen') #color='firebrick',+0.27#in m #ESD t  
		#axs[1].plot(secs/3600,esd1*1e6) #+0.27#in m #ESD t
		#axs[1].plot(secs[ymin]/3600, esd[ymin]*1e6, "o", label="min", color='r') #plot min
		#axs[1].plot(secs[ymax]/3600, esd[ymax]*1e6, "o", label="max", color='b') #plot max
#		#axs[2].plot(secs/3600,G,label=i.name) #abs(uu)*1.414  u change with time alpha=0.3,

		axs[1].plot(secs/3600,-breakup,label=i.name)
#		axs[2].plot(secs/3600,uu,label=i.name) #,color='seagreen') #color='firebrick',label="model"    u change with time alpha=0.3,  #abs(uu)*sqrt(2)
		#axs[3].plot(secs/3600,G,label=i.name) #secs/3600 G  /(uu+0.05)

		axs[2].plot(secs/3600,sinking,label=i.name)

		#axs[2].plot(Gsize,esd*1e6,alpha=0.3)  #totlpmsize/sqrt(Gsize) 	#ESD vs G
	
		#axs[3].scatter(ctot,esd*1e6, alpha=0.3) # 	G,totlpm,	#SPMC vs G

		axs[3].plot(secs/3600,dD,label=i.name) #xsize   plotting compound of all dD terms
		#axs[3].plot(secs/3600,coagulation) #,label=i.name
		#axs[3].plot(secs/3600,-breakup) #,label=i.name
		#axs[3].plot(secs/3600,sinking,label=i.name) #,label=i.name
		axs[4].plot(secs/3600,diffset,label=i.name) #xsize   plotting compound of all dD terms
		axs[5].plot(secs/3600,ws*1000,label=i.name) #xsize   plotting compound of all dD terms

		#axs[3].plot(secs/3600,Dbreakup,label=i.name)
		#axs[4].plot(secs/3600,Dsinking,label=i.name)
		#axs[4].plot(secs/3600,Derosion,label=i.name)#mu=1.1e-3	#Ws vs t
		
		ptime=secs[tslice]/3600.
		array = ptime
		
##get index for value/values to compare with spmc data/size data
		for value in tsmdatat:#tmzspmjuly[:,0]+7: #winterspmc[:,0]:#datatime: #s1datatimecorr: #[0:84]: #time in spmc data
			value = find_nearest (array, value)
			#print(find_nearest(array, value))
			Gidx.append(value) 	#append theoretical G index at certain t
			#Gidx=np.array(Gidx)
		for values in d50datat: #tmzd50july[:,0]+7:#winteresd[:,0]:#s1datatime: 	#time in size data
			values=find_nearest (array, values)
			Gsidx.append(values)	 
			#Gsidx=np.array(Gsidx)
		for valuef in flowdatat: #time in flow data
			valuef=find_nearest (array, valuef)
			Gsidf.append(valuef)	
		totlpmi=[]
		esdi=[]
		flowi=[]
		flowall=[]
		mse=[]
		mses=[]
		msef=[]
		errorsum=[]
		msei=[]
		msesi=[]
		test=[]
		test1=[]
		for i2 in Gidx:
			a=ctot[i2] #changed to ctot #totlpm[i2] totlpm at certain time for comparing #changed to 0 not using
			totlpmi.append(a)
		
		for j in Gsidx:
			c=esd[j] #size at certain time for comparing #change to 0 not using 
			esdi.append(c)
		
			#calculate MSE

		for k in Gsidf:  #find the flow velocity
			f=uu[k] #flow velocity at certain time for comparing #change to 0 not using 
			flowi.append(f) #abs(f)*sqrt(2)  #absolute value
		sumb=0
		sumd=0
		sumf=0
			
		sumbb=0
		sumdd=0
		sumff=0

		sumbt=0
		sumdt=0
		sumft=0
		for i1 in length:
			if i1<1: #8:
				b=0
				bb=0

				#b=(dataspmc[i1]-totlpmi[i1])**2 #here #MSE absolute change
			else:
				b=sqrt(((tsmdata[i1]-totlpmi[i1])/tsmdata[i1])**2)  #comment out for not using mse
				bb=(tsmdata[i1]-totlpmi[i1])**2

				#b=sqrt(((tmzspmjuly[i1,1]-totlpmi[i1])/tmzspmjuly[i1,1])**2)
				#b=sqrt(((winterspmc[i1,1]-totlpmi[i1])/winterspmc[i1,1])**2) #relative change Fettweis
				#b=sqrt(((dataspmc[i1]-totlpmi[i1])/dataspmc[i1])**2) #relative change
			sumb=sumb+b
			sumbb=sumbb+bb

			#test.append(b)
			#test01=sum(test)/len(test)
		mse=sumb/len(length)	
		msebb=sqrt(sumbb)/len(length)	
##adding a part which 1 data point difference is tolerated.####
		nt=2
		for i1t in length[nt:-nt]:
			#if i1t<1:
			#	bt=0
			#else:
			#bt=min([(tsmdata[i1t]-totlpmi[i1t-1])**2,(tsmdata[i1t]-totlpmi[i1t])**2,(tsmdata[i1t]-totlpmi[i1t+1])**2])	
			bt=min([(tsmdata[i1t]-totlpmi[i1t-2])**2,(tsmdata[i1t]-totlpmi[i1t-1])**2,(tsmdata[i1t]-totlpmi[i1t])**2,(tsmdata[i1t]-totlpmi[i1t+1])**2,(tsmdata[i1t]-totlpmi[i1t+2])**2])					
			sumbt=sumbt+bt
		msebt=sqrt(sumbt)/len(length[nt:-nt])


		for j1 in lengths:
			#d=(s1datasize[j1]-esdi[j1]*1e6)**2
			if j1<1: #6:
				d=0
				dd=0
			else:
				d=sqrt(((d50data[j1]-esdi[j1]*1e6)/d50data[j1])**2)
				dd=(d50data[j1]-esdi[j1]*1e6)**2
				#d=sqrt(((winteresd[j1,1]-esdi[j1]*1e6)/winteresd[j1,1])**2)
				#d=sqrt(((s1datasize[j1]-esdi[j1]*1e6)/s1datasize[j1])**2)
			sumd=sumd+d
			sumdd=sumdd+dd
			#test1.append(d)
			#test11=sum(test1)/len(test1)
		mses=sumd/len(lengths)
		msesdd=sqrt(sumdd)/len(lengths)
			#secs[c]/3600, esd[c]*1e6
#tolerate but for size
		for j1t in lengths[nt:-nt]:
			#if j1t<1:
			#	dt=0
			#else:
			#dt=min([(d50data[j1t]-esdi[j1t-1]*1e6)**2,(d50data[j1t]-esdi[j1t]*1e6)**2,(d50data[j1t]-esdi[j1t+1]*1e6)**2])
			dt=min([(d50data[j1t]-esdi[j1t-2]*1e6)**2,(d50data[j1t]-esdi[j1t-1]*1e6)**2,(d50data[j1t]-esdi[j1t]*1e6)**2,(d50data[j1t]-esdi[j1t+1]*1e6)**2,(d50data[j1t]-esdi[j1t+2]*1e6)**2])						
			sumdt=sumdt+dt
		msesdt=sqrt(sumdt)/len(lengths[nt:-nt])

		for f1 in lengthf:
			#d=(s1datasize[j1]-esdi[j1]*1e6)**2
			if f1<1: #6:
				ff1=0
				ff=0
			else:  #abs(uu)*1.414
				#ff=sqrt(((flowdata[f1]-(abs(flowi[f1])*sqrt(2))))**2) #/flowdata[f1]  !changed the division part
				ff1=sqrt(((flowdata[f1]-flowi[f1])/flowdata[f1])**2) #  !changed the division part
				ff=(flowdata[f1]-flowi[f1])**2
			sumf=sumf+ff1
			sumff=sumff+ff
			#test1.append(d)
			#test11=sum(test1)/len(test1)
		msef=sumf/len(lengthf)
		mseff=sqrt(sumff)/len(lengthf)
			#secs[c]/3600, esd[c]*1e6
				
		#valley=(max(esd[ymax[1]],esd[ymax[2]])-esd[ymin[2]])/min(esd[ymax[1]],esd[ymax[2]])/0.76157 if len(ymin) >= 4 else 0 	
		valley=0 #(min(esd[ymax[2]],esd[ymax[3]])-esd[ymin[3]])/min(esd[ymax[2]],esd[ymax[3]])/0.43 # 0.76157 used in bale  1.6 if len(ymin) >= 5 else 0 				#0.76157 is the valley ratio in data#*1e6   #0.43 (smaller peak - middle)/ lower value 
			#(max(esd[ymax])- max(esd[ymin]))*1e6 #.max#-esd[b].max
		#errorsum=mse+mses+(1-valley)
		
		#errorsum=mse+mses+abs(1-valley)
		minutesfrom=5*60 #7*60





		#print (err_umax)#msef (len(esd), min(esd), max(esd))
		err_dmax=abs(max(d50data)-max(esd[minutesfrom:])*1e6)#/max(d50data)
		err_dmin=abs(min(d50data)-min(esd[minutesfrom:])*1e6)#/min(d50data)
		err_dmean=abs(mean(d50data)-mean(esd[minutesfrom:])*1e6)#/mean(d50data)
		err_drange=abs(((max(d50data)-min(d50data))-(max(esd[minutesfrom:])*1e6)-min(esd[minutesfrom:])*1e6))#/(max(d50data)-min(d50data))
		err_cmax=abs(max(tsmdata)-max(ctot[minutesfrom:]))#/max(tsmdata)
		err_cmin=abs(min(tsmdata)-min(ctot[minutesfrom:]))#/min(tsmdata)
		err_cmean=abs(mean(tsmdata)-mean(ctot[minutesfrom:]))#/mean(tsmdata)
		err_crange=abs(((max(tsmdata)-min(tsmdata))-(max(ctot[minutesfrom:]))-min(ctot[minutesfrom:])))#/(max(tsmdata)-min(tsmdata))
		errorsum=(abs(err_dmax)+abs(err_dmean)+abs(err_drange)+abs(err_cmax)+abs(err_cmean)+abs(err_crange))
		print (i.name)
		#print ('errorsum with calculating range is', errorsum, "mse of SPMC is", mse, "msesize is", mses)
		#print ("mse of SPMC is", msebb, "msesize is", msesdd)
		print (msebb, msebt)
		print (msesdd, msesdt)
		stdd=abs((np.std(ctot[minutesfrom:])-np.std(tsmdata))/np.std(tsmdata))
		stdc=abs((np.std(esd[minutesfrom:]*1e6)-np.std(d50data))/np.std(d50data))

		stdd1=abs(np.std(ctot[minutesfrom:])-np.std(tsmdata))
		stdc1=abs(np.std(esd[minutesfrom:])*1e6-np.std(d50data))

		print ((err_cmax+err_cmin)/mean(tsmdata),(err_dmax+err_dmin)/mean(d50data))
	#	print ('stdc is', stdc,'stdd is', stdd) 
	#	print ('abs stdc is', stdc1,'abs stdd is', stdd1)
		#a=esd[minutesfrom:]
		#A=np.array(esd[minutesfrom:])
		#maximum_indices = np.where(A==max(a))
		A=find_max(esd[minutesfrom:])/60+minutesfrom/60 #max in d50
		diffmaxd50=(max(esd[minutesfrom:])*1e6-max(d50data)) #abs
		lagmaxd50=(A-d50datat[find_max(d50data)]) #abs

		B=np.argmin(ctot[minutesfrom:])/60+minutesfrom/60 #min in spmc
		diffminc=(min(ctot[minutesfrom:])-min(tsmdata[8:])) #abs
		lagminc=(B-tsmdatat[8+np.argmin(tsmdata[8:])]) #abs

		####flow error part############################!!!#####here!!!!!####
		#err_umax=(max(d50data)-max(flowi))/max(d50data)
		#err_umin=(min(d50data)-min(flowi))/min(d50data)
###careful with new dataset########################################
#		Fmin=np.argmin(flowi)/60+minutesfrom/60 
#		diffminf=(min(flowi)-min(flowdata/100)) #abs
#		lagminf=(Fmin-flowdatat[np.argmin(flowdata)]) #abs

#		Fmax=find_max(uu)/60+12 #+minutesfrom/60 
#		diffmaxf=(max(uu)-max(flowdata)) #abs  /100 for bale data
#		lagmaxf=(Fmax-flowdatat[np.argmax(flowdata)]) #abs
#		if abs(lagmaxf)>=4:
#			lagmaxf=lagmaxf/lagmaxf*abs(lagmaxf)%4
#			if abs(lagmaxf)>=2:
#				lagmaxf=-lagmaxf/lagmaxf*(4-abs(lagmaxf))

		#print ("flow error is",msef)#(len(esd), min(esd), max(esd))

		#print (A,d50datat[find_max(d50data)],abs(A-d50datat[find_max(d50data)]))
		#print (max(esd[minutesfrom:])*1e6,max(d50data),abs(max(esd[minutesfrom:])*1e6-max(d50data)))
################if to add the lag to the error sum function?#######################
		#print (diffminc, lagminc)
		#print (diffmaxd50,lagmaxd50)
		#print (diffminf, lagminf)
		#print (diffmaxf, lagmaxf)
#####################3here!!!!!!!!!!!!!!!!#########################################
		#print ("relative error c", (abs(err_cmax)+abs(err_crange)+abs(err_cmean))/3,"relative error d", (abs(err_dmax)+abs(err_drange)+abs(err_dmean))/3,"error sum", errorsum)

		#print ("relative error SPMC=",mse, "relative error ESD=",mses, "sum=",mse+mses, "valley=",valley, "mins(4) = ",len(ymin),"errorsum=",errorsum) #"SumMSE=",mse+mses
	#print (filelist)

###writing to file about error#############
#		f=open('relativeerror.txt', 'a') #a+
		#parts=re.split('(\d+)',i.name)
#		parts=re.split('-',i.name)
#		print (parts)
#		f.write(str(filenumber)+" ")
#		f.write(sys.argv[1]+" ")
#		f.write(parts[1]+" ") #" "
#		f.write(parts[3]+" ") #" "
#		f.write(parts[5]+" ") #" "
		##f.write(jj+" ")
		##f.write(kk+" ")
		
##		f.write(str(abs(err_cmax))+" ")
##		f.write(str(abs(err_crange))+" ")
##		f.write(str(abs(err_cmean))+" ")
##		f.write(str(abs(err_dmax))+" ")
##		f.write(str(abs(err_drange))+" ")
##

#		f.write(str(abs(err_dmean))+" ")

#		f.write(str(mse)+" ")
#		f.write(str(mses)+" ")
#		f.write(str(mse+mses)+" ")
#		f.write(str(valley)+" ")
#		f.write(str(len(ymin))+" ") #4 for valley true if run for 8 h
#		f.write(str(min(esd[720:]))+" ") #120

#		minutesfromold=4*60 #720
#		f.write(str(max(esd[minutesfromold:]))+" ")
#		f.write(str(mean(esd[minutesfromold:]))+" ")
#		f.write(str(min(ctot[minutesfromold:]))+" ") #totlpm
#		f.write(str(max(ctot[minutesfromold:]))+" ")
#		f.write(str(mean(ctot[minutesfromold:]))+" ")
#		f.write(str(min(abs(uu[minutesfromold:])*sqrt(2)))+" ")
#		f.write(str(max(abs(uu[minutesfromold:])*sqrt(2)))+" ")
#		f.write(str(mean(abs(uu[minutesfromold:])*sqrt(2)))+" ")
#		f.write(str(errorsum)+'\n')#'\n'
#		f.flush()
#		f.close()
		#print("closed")
		##axs[0].scatter(datatime, tstarttotlpmi,label=sys.argv[1]+ii+sys.argv[3]+jj)   #scatter of modelled SPMC
		##axs[1].scatter(s1datatime,1e6*np.array(esdi),label=sys.argv[1]+ii+sys.argv[3]+jj)


#		with open ("result"+str(filenumber)+".dat", 'wb') as f0:
#			pickle.dump([str(filenumber),tsmdata,totlpmi],f0)

#		f=open('results.txt', 'a') #a+
##		f.write(str(filenumber)+" ")
#		f.write(str(tsmdata)+" ")
#		f.write(str(totlpmi)+'\n')
#		f.flush()
#		f.close()

		filenumber +=1
	except OSError:
		continue
##from here plotting data
#axs[0].plot(datatime[8:],dataspmc[8:],'k^',ms=5,label='total lit, data') 		#dataSPMCvs	t

plt.rcParams.update({'font.size': 20})
#axs[0].fill_between(tsmdatat,tsmdatamax,tsmdatamin,alpha=0.30,color='grey')  #linewidth=3,
#axs[0].plot(tsmdatat,tsmdata,'o',color='grey',alpha=0.5) #,label='F&L TMZ July'
axs[0].yaxis.set_tick_params(labelsize=20)
axs[1].yaxis.set_tick_params(labelsize=20)
axs[2].xaxis.set_tick_params(labelsize=20)
axs[2].yaxis.set_tick_params(labelsize=20)
#axs[3].yaxis.set_tick_params(labelsize=20)
#axs[3].xaxis.set_tick_params(labelsize=20)


axs[0].set_ylabel('coagulation', fontsize=20 ) #'SPMC (g m$^{-3}$)'
axs[0].set_ylim(0,)


#axs[0].scatter(B,min(ctot[minutesfrom:]))
#axs[0].scatter(tsmdatat[8+np.argmin(tsmdata[8:])],min(tsmdata[8:]))

#axs[1].plot(d50datat,d50data,'o',color='grey',alpha=0.5) #,label='F&L TMZ July'

#axs[1].plot(d50datat,d50datamax,'o',color='grey',alpha=0.5)
#axs[1].plot(d50datat,d50datamin,'o',color='grey',alpha=0.5) 


#axs[1].fill_between(d50datat,d50datamax,d50datamin,alpha=0.30, color='grey') #linewidth=3,

#axs[1].set_ylim(,0)

#axs[1].scatter(A,max(esd[minutesfrom:])*1e6)
#axs[1].scatter(d50datat[find_max(d50data)],max(d50data))

#axs[3].plot(d50datat,d50data/100-0.52,'o') #for compound

#axs[1].plot(s1datatime[6:],s1datasize[6:],'g^',ms=5.,label='ESD, data') 			#dataESD vs	t
#axs[1].set_xlim(0,8)

#dataindex = np.diff(np.sign(np.diff(s1datasize[6:]))).nonzero()[0] + 1+6               # local min & max
#datamin = (np.diff(np.sign(np.diff(s1datasize[6:]))) > 0).nonzero()[0] + 1+6         # local min
#datamax = (np.diff(np.sign(np.diff(s1datasize[6:]))) < 0).nonzero()[0] + 1+6         # local max# +1 due to
#datavalley=(min(s1datasize[datamax])- min(s1datasize[datamin]))
#print "datavalley=",datavalley


#axs[1].set_xlabel('time (h)')
axs[1].set_ylabel('breakup', fontsize=20 )  #size ($\mu m$)
#axs[1].set_xlim(0,16)
#axs[4].plot(s1datasize[6:],s1datasize[6:],'k')
Gmodel=[]
Dmodel=[]

#for i in Gsidx:
#	Gi=Gsize[i] #Gsize at certain time for plotting D vs G

	#Di=esd[i+17] #esd[i+17]+0.28h*60
	#Di=esd[i]
	#Dmodel.append(Di)
#	Gmodel.append(Gi)

#axs[2].scatter(Gmodel,s1datasize,label='data') #[5:84]  				#dataESD vs	G
#axs[2].plot(flowdatat,flowdata,'o',label='baled50flow') #'F&L TMZ July'  /100 #divided by 100 bale data



#axs[2].scatter(Fmax,max(uu))
#axs[2].scatter(flowdatat[np.argmax(flowdata)],max(flowdata))  #/100 

#axs[2].plot(flowdatat,flowdata,'o',color='grey',alpha=0.5,label='OSZ May') #'F&L TMZ July' /100 for labdata

#axs[2].fill_between(flowdatat,flowdatamax,flowdatamin,alpha=0.30, color='grey',label='data') #'tab:orange'  #linewidth=3,
axs[5].set_xlabel('time (h)',fontsize=20 )#('G')

axs[2].set_ylabel('sinking',fontsize=20 )#('ESD ($\mu$m)') velocity  #'flow (m s$^{-1}$)'
#axs[3].set_ylabel('shear G (s$^{-1}$)', fontsize=20 )
axs[3].set_ylabel('dD',fontsize=20 )#('ESD ($\mu$m)') velocity  #'flow (m s$^{-1}$)'
axs[4].set_ylabel('diffset',fontsize=20 )#('ESD ($\mu$m)') velocity  #'flow (m s$^{-1}$)'
axs[5].set_ylabel('w$_s$ (mm/s)',fontsize=20 )#('ESD ($\mu$m)') velocity  #'flow (m s$^{-1}$)'

plt.legend()

#axs[3].set_xlabel('SPMC (g m$^{-3}$)')
#axs[3].set_ylabel('ESD ($\mu$m)')
#axs[4].set_ylabel('dD (m s$^{-1}$)') #'ws (m s$^{-1}$)' ('modelled SPMC')# ('modelled ESD ($\mu$m)')

plt.show()
