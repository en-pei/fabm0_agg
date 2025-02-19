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

filelist=[]
path='/gpfs/home/lie/setups/bale2002/par'
#input from runpost.py to iterate through all layers
lay=sys.argv[1] #5 #16 #4 #16 #6 #16 #layer for SPM concentration
laysize=sys.argv[2]# 5 #16 #14 #9 #layer for ESD
print(lay,laysize,'lay')

h=3600./60 #/120. output was per 120s
tstart=int(0*h)#int(6.*h)
tend=tstart+int(39*h)#tstart+int(8*h)
tslice=slice(tstart,tend)

datadir = '/gpfs/home/lie/setups/bale2002/Plot/data'
data = asarray(numpy.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=', '))
sizedata = asarray(numpy.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))
s1datatime=sizedata[:,0]-5#-0.25
s1datasize=sizedata[:,1] #data is in micrometer
datatime=data[:,0]-5.#-0.25
dataspmc=data[:,2]

winteresd=asarray(numpy.loadtxt('%s/esd_feb_data.txt'%datadir,delimiter=', '))
summeresd=asarray(numpy.loadtxt('%s/esd_july_data.txt'%datadir,delimiter=', '))
winterspmc=asarray(numpy.loadtxt('%s/spmc_feb_data.txt'%datadir,delimiter=', '))
summerspmc=asarray(numpy.loadtxt('%s/spmc_july_data.txt'%datadir,delimiter=', '))

winteresdm=asarray(numpy.loadtxt('%s/esd_feb_model.txt'%datadir,delimiter=', '))
summeresdm=asarray(numpy.loadtxt('%s/esd_july_model.txt'%datadir,delimiter=', '))
winterspmcm=asarray(numpy.loadtxt('%s/spmc_feb_model.txt'%datadir,delimiter=', '))
summerspmcm=asarray(numpy.loadtxt('%s/spmc_july_model.txt'%datadir,delimiter=', '))

bloomspmc=asarray(numpy.loadtxt('%s/bloomspm.txt'%datadir,delimiter=', '))
regularspmc=asarray(numpy.loadtxt('%s/regularspm.txt'%datadir,delimiter=', '))
bloomd50=asarray(numpy.loadtxt('%s/bloomd50.txt'%datadir,delimiter=', '))
regulard50=asarray(numpy.loadtxt('%s/regulard50.txt'%datadir,delimiter=', '))

tmzspmjuly=asarray(numpy.loadtxt('%s/tmzspmjuly.txt'%datadir,delimiter=', '))
tmzd50july=asarray(numpy.loadtxt('%s/tmzd50july.txt'%datadir,delimiter=', '))
tmzflowjuly=asarray(numpy.loadtxt('%s/tmzflowjuly.txt'%datadir,delimiter=', '))
tmzspmfeb=asarray(numpy.loadtxt('%s/tmzspmfeb.txt'%datadir,delimiter=', '))
tmzd50feb=asarray(numpy.loadtxt('%s/tmzd50feb.txt'%datadir,delimiter=', '))
tmzflowfeb=asarray(numpy.loadtxt('%s/tmzflowfeb.txt'%datadir,delimiter=', '))

oszspmmay=asarray(numpy.loadtxt('%s/oszspmmay.txt'%datadir,delimiter=', '))
oszd50may=asarray(numpy.loadtxt('%s/oszd50may.txt'%datadir,delimiter=', '))
oszflowmay=asarray(numpy.loadtxt('%s/oszflowmay.txt'%datadir,delimiter=', '))
oszspmfeb=asarray(numpy.loadtxt('%s/oszspmfeb.txt'%datadir,delimiter=', '))
oszd50feb=asarray(numpy.loadtxt('%s/oszd50feb.txt'%datadir,delimiter=', '))
oszflowfeb=asarray(numpy.loadtxt('%s/oszflowfeb.txt'%datadir,delimiter=', '))

bale3500spmc=asarray(numpy.loadtxt('%s/spmc_bale3500.txt'%datadir,delimiter=', '))
bale3500spmcf=asarray(numpy.loadtxt('%s/flow_balespmc3500.txt'%datadir,delimiter=', '))
bale3500d50=asarray(numpy.loadtxt('%s/d50_bale3500.txt'%datadir,delimiter=', '))
bale3500d50f=asarray(numpy.loadtxt('%s/flow_baled503500.txt'%datadir,delimiter=', '))


##define dataset
tsmdatat=bale3500spmc[:,0]-0.95-0.1-0.42 #datatime #tmzspmjuly[:,0]+7.5 #oszspmfeb[:,0] #tmzspmfeb[:-1,0]+4+6 #
d50datat=bale3500d50[:,0]-0.95-0.1-0.42 #s1datatime #tmzd50july[:,0]+7.5 #oszd50feb[:,0] #tmzd50feb[:,0]+4+6 #
flowdatat=bale3500d50f[:,0]-0.95-0.1-0.42 #tmzflowjuly[:,0]+7.5 #oszflowfeb[:,0] #tmzflowfeb[:,0]+4+6 #
flowdatat1=bale3500spmcf[:,0]-0.95-0.1-0.42 

tsmdata=bale3500spmc[:,1] #dataspmc #tmzspmjuly[:,1] #oszspmfeb[:,1] #tmzspmfeb[:-1,1] #
d50data=bale3500d50[:,1] #s1datasize #tmzd50july[:,1] #oszd50feb[:,1]  #tmzd50feb[:,1] #
flowdata=bale3500d50f[:,1] #tmzflowjuly[:,1] #oszflowfeb[:,1] #tmzflowfeb[:,1] #
flowdata1=bale3500spmcf[:,1] 

s1datatimecorr=list(i for i in s1datatime) #(i-0.27) #time in datasize, 2.7 h earlier
#print len(datatime) #79
length=range(len(tsmdatat)) #datatime	len(winterspmc[:,0])			#idexing spmc data
lengths=range(len(d50datat)) 	#s1datasize len(winteresd[:,0])			#indexing size data
lengthf=range(len(flowdatat)) 	#s1datasize len(winteresd[:,0])			#indexing size data


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx#, array[idx]

Gidx=[]
Gsidx=[]
Gsidf=[]
filenumber=0
#plt.rcParams.update({'font.size': 14})
fig, axs = plt.subplots(4, sharex=True)  #(5, sharex=False)
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
		phy=squeeze(ncv['npzd_phy'][:,lay])
		det=squeeze(ncv['npzd_det'][:,lay])
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
		dD=squeeze(ncv['agg_resuspensionlpm'][:,laysize])
		xsize=squeeze(ncv['agg_Xsize'][:,laysize])
#		Dbreakup=squeeze(ncv['agg_Dbreakup'][:,laysize])
#		Dcoagulation=squeeze(ncv['agg_Dcoagulation'][:,laysize])
#		Dsinking=squeeze(ncv['agg_sinkinglpm'][:,laysize])
#		Derosion=squeeze(ncv['agg_dflux'][:])
		uu = squeeze(ncv['u'][:,laysize]) #flow velocity
		#ws2 = squeeze(ncv['agg_sinkinglpm'][:,lay]) #agg_Breakup
		#ws = -ws1 +ws2
		aindex = np.diff(np.sign(np.diff(esd))).nonzero()[0] + 1             # local min & max
		ymin = (np.diff(np.sign(np.diff(esd))) > 0).nonzero()[0] + 1         # local min
		ymax = (np.diff(np.sign(np.diff(esd))) < 0).nonzero()[0] + 1         # local max# +1 
		axs[0].plot(secs/3600,ctot,label=i.name) #totlpm    #+0.27 #SPMC	vs t
		#axs[0].plot(secs/3600,phy+det,label=i.name) #totlpm    #+0.27 #SPMC	vs t
		#axs[0].plot(secs/3600,ctot1) #totlpm    #+0.27 #SPMC	vs t
		axs[1].plot(secs/3600,esd*1e6,label=i.name) #+0.27#in m #ESD t
		#axs[1].plot(secs/3600,esd1*1e6) #+0.27#in m #ESD t
		#axs[1].plot(secs[ymin]/3600, esd[ymin]*1e6, "o", label="min", color='r') #plot min
		#axs[1].plot(secs[ymax]/3600, esd[ymax]*1e6, "o", label="max", color='b') #plot max
		axs[2].plot(secs/3600,G,label=i.name) #abs(uu)*1.414  u change with time alpha=0.3,
		#axs[2].plot(Gsize,esd*1e6,alpha=0.3)  #totlpmsize/sqrt(Gsize) 	#ESD vs G
	
		#axs[3].scatter(ctot,esd*1e6, alpha=0.3) # 	G,totlpm,	#SPMC vs G

		axs[3].plot(secs/3600,dD,label=i.name) #xsize   plotting compound
		#axs[3].plot(secs/3600,coagulation,label=i.name) 
		#axs[3].plot(secs/3600,-breakup,label=i.name) 
		#axs[3].plot(secs/3600,sinking,label=i.name) 


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
			flowi.append(f)
		sumb=0
		sumd=0
		sumf=0
			
		for i1 in length:
			if i1<8:
				b=0
				#b=(dataspmc[i1]-totlpmi[i1])**2 #here #MSE absolute change
			else:
				b=sqrt(((tsmdata[i1]-totlpmi[i1])/tsmdata[i1])**2)  #comment out for not using mse
				#b=sqrt(((tmzspmjuly[i1,1]-totlpmi[i1])/tmzspmjuly[i1,1])**2)
				#b=sqrt(((winterspmc[i1,1]-totlpmi[i1])/winterspmc[i1,1])**2) #relative change Fettweis
				#b=sqrt(((dataspmc[i1]-totlpmi[i1])/dataspmc[i1])**2) #relative change
			sumb=sumb+b
			#test.append(b)
			#test01=sum(test)/len(test)
			mse=sumb/len(length)				
		for j1 in lengths:
			#d=(s1datasize[j1]-esdi[j1]*1e6)**2
			if j1<6:
				d=0
			else:
				d=sqrt(((d50data[j1]-esdi[j1]*1e6)/d50data[j1])**2)
				#d=sqrt(((winteresd[j1,1]-esdi[j1]*1e6)/winteresd[j1,1])**2)
				#d=sqrt(((s1datasize[j1]-esdi[j1]*1e6)/s1datasize[j1])**2)
			sumd=sumd+d
			#test1.append(d)
			#test11=sum(test1)/len(test1)
			mses=sumd/len(lengths)
			#secs[c]/3600, esd[c]*1e6

		for f1 in lengthf:
			#d=(s1datasize[j1]-esdi[j1]*1e6)**2
			if f1<1: #6:
				ff=0
			else:
				ff=sqrt(((flowdata[f1]-(abs(flowi[f1])*sqrt(2)))/flowdata[f1])**2)
			sumf=sumf+ff
			#test1.append(d)
			#test11=sum(test1)/len(test1)
			msef=sumf/len(lengthf)
			#secs[c]/3600, esd[c]*1e6
				
		#valley=(max(esd[ymax[1]],esd[ymax[2]])-esd[ymin[2]])/min(esd[ymax[1]],esd[ymax[2]])/0.76157 if len(ymin) >= 4 else 0 	
		valley=(min(esd[ymax[2]],esd[ymax[3]])-esd[ymin[3]])/min(esd[ymax[2]],esd[ymax[3]])/1.6 if len(ymin) >= 5 else 0 				#0.76157 is the valley ratio in data#*1e6
			#(max(esd[ymax])- max(esd[ymin]))*1e6 #.max#-esd[b].max
		#errorsum=mse+mses+(1-valley)
		
		#errorsum=mse+mses+abs(1-valley)
		minutesfrom=4*60 #7*60
		print ("flow error is",msef)#(len(esd), min(esd), max(esd))
		err_dmax=(max(d50data)-max(esd[minutesfrom:])*1e6)/max(d50data)
		err_dmin=(min(d50data)-min(esd[minutesfrom:])*1e6)/min(d50data)
		err_dmean=(mean(d50data)-mean(esd[minutesfrom:])*1e6)/mean(d50data)
		err_drange=(((max(d50data)-min(d50data))-(max(esd[minutesfrom:])*1e6)-min(esd[minutesfrom:])*1e6))/(max(d50data)-min(d50data))
		err_cmax=(max(tsmdata)-max(ctot[minutesfrom:]))/max(tsmdata)
		err_cmin=(min(tsmdata)-min(ctot[minutesfrom:]))/min(tsmdata)
		err_cmean=0 #!!!!!!!#(mean(tsmdata)-mean(ctot[minutesfrom:]))/mean(tsmdata)
		err_crange=(((max(tsmdata)-min(tsmdata))-(max(ctot[minutesfrom:]))-min(ctot[minutesfrom:])))/(max(tsmdata)-min(tsmdata))
		errorsum=(abs(err_dmax)+abs(err_dmean)+abs(err_drange)+abs(err_cmax)+abs(err_cmean)+abs(err_crange))
		print ("relative error c", (abs(err_cmax)+abs(err_crange)+abs(err_cmean))/3,"relative error d", (abs(err_dmax)+abs(err_drange)+abs(err_dmean))/3,"error sum", errorsum)

		#print ("relative error SPMC=",mse, "relative error ESD=",mses, "sum=",mse+mses, "valley=",valley, "mins(4) = ",len(ymin),"errorsum=",errorsum) #"SumMSE=",mse+mses
	#print (filelist)
		f=open('relativeerror.txt', 'a') #a+
		#parts=re.split('(\d+)',i.name)
		parts=re.split('-',i.name)
		print (parts)
		f.write(str(filenumber)+" ")
		f.write(sys.argv[1]+" ")
		f.write(parts[1]+" ") #" "
		f.write(parts[3]+" ") #" "
		f.write(parts[5]+" ") #" "
		##f.write(jj+" ")
		##f.write(kk+" ")
		
#		f.write(str(abs(err_cmax))+" ")
#		f.write(str(abs(err_crange))+" ")
#		f.write(str(abs(err_cmean))+" ")
#		f.write(str(abs(err_dmax))+" ")
#		f.write(str(abs(err_drange))+" ")
#		f.write(str(abs(err_dmean))+" ")

		f.write(str(mse)+" ")
		f.write(str(mses)+" ")
		f.write(str(mse+mses)+" ")
		f.write(str(valley)+" ")
		f.write(str(len(ymin))+" ") #4 for valley true if run for 8 h
		f.write(str(min(esd[720:]))+" ") #120

		minutesfromold=4*60 #720
		f.write(str(max(esd[minutesfromold:]))+" ")
		f.write(str(mean(esd[minutesfromold:]))+" ")
		f.write(str(min(ctot[minutesfromold:]))+" ") #totlpm
		f.write(str(max(ctot[minutesfromold:]))+" ")
#		f.write(str(mean(ctot[minutesfromold:]))+" ")
		f.write(str(min(abs(uu[minutesfromold:])*sqrt(2)))+" ")
		f.write(str(max(abs(uu[minutesfromold:])*sqrt(2)))+" ")
		f.write(str(mean(abs(uu[minutesfromold:])*sqrt(2)))+" ")
		f.write(str(errorsum)+'\n')#'\n'
		f.flush()
		f.close()
		#print("closed")
		##axs[0].scatter(datatime, tstarttotlpmi,label=sys.argv[1]+ii+sys.argv[3]+jj)   #scatter of modelled SPMC
		##axs[1].scatter(s1datatime,1e6*np.array(esdi),label=sys.argv[1]+ii+sys.argv[3]+jj)


		with open ("result"+str(filenumber)+".dat", 'wb') as f0:
			pickle.dump([str(filenumber),tsmdata,totlpmi],f0)

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

axs[0].plot(tsmdatat,tsmdata,'o') #,label='F&L TMZ July'
axs[0].set_ylabel('SPMC (g m$^{-3}$)')
#axs[0].set_xlim(0,8)

axs[1].plot(d50datat,d50data,'o') #,label='F&L TMZ July'
#axs[3].plot(d50datat,d50data/100-0.52,'o') #for compound

#axs[1].plot(s1datatime[6:],s1datasize[6:],'g^',ms=5.,label='ESD, data') 			#dataESD vs	t
#axs[1].set_xlim(0,8)

#dataindex = np.diff(np.sign(np.diff(s1datasize[6:]))).nonzero()[0] + 1+6               # local min & max
#datamin = (np.diff(np.sign(np.diff(s1datasize[6:]))) > 0).nonzero()[0] + 1+6         # local min
#datamax = (np.diff(np.sign(np.diff(s1datasize[6:]))) < 0).nonzero()[0] + 1+6         # local max# +1 due to
#datavalley=(min(s1datasize[datamax])- min(s1datasize[datamin]))
#print "datavalley=",datavalley


#axs[1].set_xlabel('time (h)')
axs[1].set_ylabel('mean ESD ($\mu m$)')
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

axs[2].plot(flowdatat,flowdata/100,'o',label='baled50flow') #'F&L TMZ July'
axs[2].plot(flowdatat1,flowdata1/100,'o',label='balespmcflow') #'F&L TMZ July'
axs[2].set_xlabel('time (h)')#('G')
axs[2].set_ylabel('flow velocity (m s$^{-1}$)')#('ESD ($\mu$m)')
plt.legend()

#axs[3].set_xlabel('SPMC (g m$^{-3}$)')
#axs[3].set_ylabel('ESD ($\mu$m)')
#axs[4].set_ylabel('dD (m s$^{-1}$)') #'ws (m s$^{-1}$)' ('modelled SPMC')# ('modelled ESD ($\mu$m)')

plt.show()
