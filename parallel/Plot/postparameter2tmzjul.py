from matplotlib import rc
from matplotlib import ticker, cm
import netCDF4
from pylab import *
from scipy.interpolate import interp1d
import cftime #netcdftime
import numpy	
import os
import sys
import subprocess

###for slurm_agg1.sh!!!!######
####code for FETTWEIS field range identification#####
##this version not relative error but absolute eror, avoid division#######
##toleration to values around 2 data point plus and minus#####


#from glob import glob
####to delete file while generating
filelist=[]
#path='/home/enpei/mossco/setups/bale2002/paramstudy'
path='/gpfs/home/lie/setups/bale2002/testtmzjul' #test1'

lay=6 #2 #16 #6 #4 #16 #layer for SPM concentration
laysize=6 #2 #9 #6 #4 #9 #layer for ESD


h=3600./60 #/120. output was per 120s
tstart=int(0*h)#int(6.*h)
tend=tstart+int(8*h) #tstart+int(39*h)#
tslice=slice(tstart,tend)

#datadir = '/home/enpei/mossco/setups/bale2002/Plot/data'
datadir='/gpfs/home/lie/setups/bale2002/Plot/data'

###LAB DATA#######
#data = asarray(np.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=', '))
#sizedata = asarray(np.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))
#s1datatime=sizedata[:,0]-5
#s1datasize=sizedata[:,1] #data is in micrometer
#datatime=data[:,0]-5.
#dataspmc=data[:,2]
###LAB DATA#######
#bale3500spmc=asarray(numpy.loadtxt('%s/spmc_bale3500.txt'%datadir,delimiter=', '))
##bale3500spmcf=asarray(numpy.loadtxt('%s/flow_balespmc3500.txt'%datadir,delimiter=', '))
#bale3500d50=asarray(numpy.loadtxt('%s/d50_bale3500.txt'%datadir,delimiter=', '))
#bale3500d50f=asarray(numpy.loadtxt('%s/flow_baled503500.txt'%datadir,delimiter=', '))

###dataset for oszmay############
#oszspmmay=asarray(numpy.loadtxt('%s/oszspmmay.txt'%datadir,delimiter=', '))
#oszd50may=asarray(numpy.loadtxt('%s/oszd50may.txt'%datadir,delimiter=', '))
#oszflowmay=asarray(numpy.loadtxt('%s/oszflowmay.txt'%datadir,delimiter=', '))

tmzspmjuly=asarray(numpy.loadtxt('%s/tmzspmjuly.txt'%datadir,delimiter=', '))
tmzd50july=asarray(numpy.loadtxt('%s/tmzd50july.txt'%datadir,delimiter=', '))
tmzflowjuly=asarray(numpy.loadtxt('%s/tmzflowjuly.txt'%datadir,delimiter=', '))

##define dataset
####tsmdatat=bale3500spmc[:,0]-0.95-0.1 #tmzspmjuly[:,0]+7.5 #tmzspmfeb[:-1,0]+4+6+0.5+8.5 #oszspmmay[:,0]+7 #oszspmfeb[:,0]+7 # #datatime #####
########################now oszmay########################################3

tsmdatat=tmzspmjuly[:,0]+7.5 ##bale3500spmc[8:,0]-0.95-0.1 # oszspmmay[:,0]+6-0.3 #+7 #-0.42 #datatime #oszspmfeb[:,0] #tmzspmfeb[:-1,0]+4+6 #
d50datat=tmzd50july[:,0]+7.5 #bale3500d50[8:,0]-0.95-0.1  # oszd50may[:,0]+6-0.3 #+7 #-0.42 #s1datatime #oszd50feb[:,0] #tmzd50feb[:,0]+4+6 #
flowdatat=tmzflowjuly[:,0]+7.5 #bale3500d50f[8:,0]-0.95-0.1 # oszflowmay[:,0]+6-0.3 #+7 #-0.42 #oszflowfeb[:,0] #tmzflowfeb[:,0]+4+6 #
#flowdatat1=bale3500spmcf[:,0]-0.95-0.1#-0.42 

tsmdata=tmzspmjuly[:,1] #bale3500spmc[8:,1] #oszspmmay[:,1] #dataspmc #oszspmfeb[:,1] #tmzspmfeb[:-1,1] #
d50data=tmzd50july[:,1] #bale3500d50[8:,1] #oszd50may[:,1] #s1datasize #oszd50feb[:,1]  #tmzd50feb[:,1] #
flowdata=tmzflowjuly[:,1] #bale3500d50f[8:,1] # oszflowmay[:,1] #oszflowfeb[:,1] #tmzflowfeb[:,1] #
#flowdata1=bale3500spmcf[:,1] 

#print (tsmdatat)
#exit() 

#s1datatimecorr=list(i for i in s1datatime) #(i-0.27) #time in datasize, 2.7 h earlier
#print len(datatime) #79
length=range(len(tsmdatat))  #range(len(datatime)) 				#idexing spmc data
lengths=range(len(d50datat)) #range(len(s1datasize)) 				#indexing size data

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx#, array[idx]

Gidx=[]
Gsidx=[]

pathsub=sys.argv[1]
ncfile='bale2002.nc'
nc=netCDF4.Dataset(ncfile)
ncv=nc.variables
utime=cftime.utime(ncv['time'].units) #netcdftime.utime(ncv['time'].units)
time=utime.num2date(ncv['time'][:])
secs=ncv['time'][:]
if len(secs)<481:
	exit() #skip run time smaller ones
#for i in os.scandir(path):
#	i.is_dir()
#	filelist.append(i.name)
#	try:
#		ncfile=i.path+'/bale2002.nc'
#		nc=netCDF4.Dataset(ncfile)
#		ncv=nc.variables
#		utime=cftime.utime(ncv['time'].units) #netcdftime.utime(ncv['time'].units)
#		time=utime.num2date(ncv['time'][:])
#		secs=ncv['time'][:]
#		#print (len(secs))
#		if len(secs)<481:
#			continue #skip run time smaller ones
lpm = squeeze(ncv['spm_spm'][:,lay])  #[:,laysize]
agglpm = squeeze(ncv['agg_agglpm'][:,lay]) #[:,laysize]
totlpm=lpm+agglpm
lpmsize = squeeze(ncv['spm_spm'][:,laysize])  #in size layer
agglpmsize = squeeze(ncv['agg_agglpm'][:,laysize])

dsize=squeeze(ncv['agg_Dsize'][:,laysize]) #dynamical size
esd=dsize #squeeze(ncv['agg_esd'][:,laysize])  #dsize #dsize is for dynamical
###doesn't use G ####
#Gsize = squeeze(ncv['agg_G'][:,laysize])
#G = squeeze(ncv['agg_G'][:,lay])
	#	ws = squeeze(ncv['agg_resuspensionlpm'][:,lay]) #(ncv['agg_ws'][:,lay])
		#ws2 = squeeze(ncv['agg_sinkinglpm'][:,lay]) #agg_Breakup
		#ws = -ws1 +ws2
aindex = np.diff(np.sign(np.diff(esd))).nonzero()[0] + 1             # local min & max
ymin = (np.diff(np.sign(np.diff(esd))) > 0).nonzero()[0] + 1         # local min
ymax = (np.diff(np.sign(np.diff(esd))) < 0).nonzero()[0] + 1         # local max# +1 
ptime=secs[tslice]/3600.
array = ptime
					#get index for value/values to compare with spmc data/size data
for value in tsmdatat: #datatime: #s1datatimecorr: #[0:84]: #time in spmc data
	value = find_nearest (array, value)
			#print(find_nearest(array, value))
	Gidx.append(value) 	#append theoretical G index at certain t
			#Gidx=np.array(Gidx)
for values in d50datat: #s1datatime: 	#time in size data
	values=find_nearest (array, values)
	Gsidx.append(values)	 
			#Gsidx=np.array(Gsidx)
totlpmi=[]
esdi=[]
mse=[]
mses=[]
errorsum=[]
msei=[]
msesi=[]
test=[]
test1=[]
for i2 in Gidx:
	a=totlpm[i2] #totlpm at certain time for comparing
	totlpmi.append(a)
			#totlpmi=np.array(totlpmi)
	
for j in Gsidx:
	c=esd[j] #size at certain time for comparing
	esdi.append(c)
			#esdi=np.array(esdi)
		
			#calculate MSE
sumb=0
sumd=0
sumbb=0
sumdd=0
sumbt=0
sumdt=0
#print (length)
###[0:25] only for oszmay###################
#######should be nothing for lab case#######
for i1 in length: #[0:25]:
	#if i1<1: #8: #8 for lab cases
	#	b=0
	#else:
	#	b=sqrt(((tsmdata[i1]-totlpmi[i1])/tsmdata[i1])**2)  #newly datathief bale
	#sumb=sumb+b
	bb=(tsmdata[i1]-totlpmi[i1])**2
	sumbb=sumbb+bb
msebb=sqrt(sumbb)/len(length)	
#mse=sumb/len(length) #[0:25]
	#print (sumb)
for j1 in lengths: #[0:25]:
	dd=(d50data[j1]-esdi[j1]*1e6)**2
	sumdd=sumdd+dd
msesdd=sqrt(sumdd)/len(lengths)
	#if j1<1: #6: #6 for lab cases
	#	d=0
	#else:
	#	d=sqrt(((d50data[j1]-esdi[j1]*1e6)/d50data[j1])**2) #newly datathief bale
	#sumd=sumd+d
#mses=sumd/len(lengths[0:25])
			#secs[c]/3600, esd[c]*1e6


nt=2
for i1t in length[nt:-nt]:
	bt=min([(tsmdata[i1t]-totlpmi[i1t-2])**2,(tsmdata[i1t]-totlpmi[i1t-1])**2,(tsmdata[i1t]-totlpmi[i1t])**2,(tsmdata[i1t]-totlpmi[i1t+1])**2,(tsmdata[i1t]-totlpmi[i1t+2])**2])					
	sumbt=sumbt+bt
msebt=sqrt(sumbt)/len(length[nt:-nt])

for j1t in lengths[nt:-nt]:
	dt=min([(d50data[j1t]-esdi[j1t-2]*1e6)**2,(d50data[j1t]-esdi[j1t-1]*1e6)**2,(d50data[j1t]-esdi[j1t]*1e6)**2,(d50data[j1t]-esdi[j1t+1]*1e6)**2,(d50data[j1t]-esdi[j1t+2]*1e6)**2])						
	sumdt=sumdt+dt
msesdt=sqrt(sumdt)/len(lengths[nt:-nt])


###valley not used for field case				
#valley=(max(esd[ymax[1]],esd[ymax[2]])-esd[ymin[2]])/min(esd[ymax[1]],esd[ymax[2]])/0.43 if len(ymin) >= 4 else 0 			#0.76157 is the valley ratio in data#*1e6
			#(max(esd[ymax])- max(esd[ymin]))*1e6 #.max#-esd[b].max
valley=0
#valley=abs((max(esd[ymax[1]],esd[ymax[2]])-esd[ymin[2]])*1e6-min(esd[ymax[1]],esd[ymax[2]])*1e6) if len(ymin) >= 4 else 0 


minutesfrom=4*60
ctot=totlpm
#from here not normalize error to avoid underestimation of max values#######
err_dmax=abs(max(d50data)-max(esd[minutesfrom:])*1e6)#/max(d50data) #normalized error from max size
err_dmin=abs(min(d50data)-min(esd[minutesfrom:])*1e6)#/min(d50data)
err_dmean=abs(mean(d50data)-mean(esd[minutesfrom:])*1e6)#/mean(d50data)
err_drange=abs(((max(d50data)-min(d50data))-(max(esd[minutesfrom:])*1e6)-min(esd[minutesfrom:])*1e6))#/(max(d50data)-min(d50data))
err_cmax=abs(max(tsmdata)-max(ctot[minutesfrom:]))#/max(tsmdata)
err_cmin=abs(min(tsmdata)-min(ctot[minutesfrom:]))#/min(tsmdata)
err_cmean=abs(mean(tsmdata)-mean(ctot[minutesfrom:]))#/mean(tsmdata)
err_crange=abs((max(tsmdata)-min(tsmdata))-(max(ctot[minutesfrom:])-min(ctot[minutesfrom:])))#/(max(tsmdata)-min(tsmdata))
rangesum=(abs(err_dmax)+abs(err_dmean)+abs(err_drange)+abs(err_cmax)+abs(err_cmean)+abs(err_crange))
stdd=abs(np.std(ctot[minutesfrom:])-np.std(tsmdata))#/np.std(tsmdata))
stdc=abs(np.std(esd[minutesfrom:]*1e6)-np.std(d50data))#/np.std(d50data))
errorsum=msebb+msesdd#+abs(valley-25)

#stdd=abs((statistics.stdev(ctot[minutesfrom:])-statistics.stdev(tsmdata))/statistics.stdev(tsmdata))
#stdc=abs((statistics.stdev(esd[minutesfrom:]*1e6)-statistics.stdev(d50data))/statistics.stdev(d50data))
 #abs(1-valley)  #for lab
#errorsum=msebb+msesdd+msebt+msesdt

f=open('../relativeerror.txt', 'a+')
#		parts=re.split('(\d+)',i.name) 
#		f.write(parts[1]+" ") #" "
#		f.write(parts[3]+" ") #" "
#		f.write(parts[5]+" ") #" "
		#f.write(jj+" ")
		#f.write(kk+" ")
f.write(pathsub+" ") #"ID"+
f.write(str(msebb)+" ")##mse  #"spmc"+  
f.write(str(msesdd)+" ") #mses ##"size"
f.write(str(msebt)+" ") #"spmc+size"
f.write(str(msesdt)+" ") #"valley"+
f.write(str(len(ymin))+" ") #4 for valley true if run for 8 h "peaknr"
f.write(str(err_dmax)+" ")
f.write(str(err_dmin)+" ")
f.write(str(err_dmean)+" ")
f.write(str(err_drange)+" ")
f.write(str(err_cmax)+" ")
f.write(str(err_cmin)+" ")
f.write(str(err_cmean)+" ")
f.write(str(err_crange)+" ")
f.write(str(rangesum)+" ")
f.write(str(stdd)+" ")
f.write(str(stdc)+" ")
f.write(str(errorsum)+'\n')#'\n'"errortot"
f.flush()
f.close()
		#axs[0].scatter(datatime, totlpmi,label=sys.argv[1]+ii+sys.argv[3]+jj)   #scatter of modelled SPMC
		#axs[1].scatter(s1datatime,1e6*np.array(esdi),label=sys.argv[1]+ii+sys.argv[3]+jj)
		#print ("relative error SPMC=",mse, "relative error ESD=",mses, "sum=",mse+mses, "valley=",valley, "mins(4) = ",len(ymin),"errorsum=",errorsum) #"SumMSE=",mse+mses
	#print (filelist)
