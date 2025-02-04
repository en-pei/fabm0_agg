#module load applications/python/3.8 
#import matplotlib
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
#from glob import glob
####to delete file while generating
filelist=[]
#path='/home/enpei/mossco/setups/bale2002/paramstudy'
path='/gpfs/home/lie/setups/bale2002/test'

lay= 16 #layer for SPM concentration
laysize= 9 #layer for ESD


h=3600./60 #/120. output was per 120s
tstart=int(0*h)#int(6.*h)
tend=tstart+int(8*h) #tstart+int(39*h)#
tslice=slice(tstart,tend)

#datadir = '/home/enpei/mossco/setups/bale2002/Plot/data'
datadir='/gpfs/home/lie/setups/bale2002/Plot/data'

data = asarray(np.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=', '))
sizedata = asarray(np.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))
s1datatime=sizedata[:,0]-5
s1datasize=sizedata[:,1] #data is in micrometer
datatime=data[:,0]-5.
dataspmc=data[:,2]

s1datatimecorr=list(i for i in s1datatime) #(i-0.27) #time in datasize, 2.7 h earlier
#print len(datatime) #79
length=range(len(datatime)) 				#idexing spmc data
lengths=range(len(s1datasize)) 				#indexing size data

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
totlpmsize=lpmsize+agglpmsize
dsize=squeeze(ncv['agg_Dsize'][:,laysize]) #dynamical size
esd=dsize #squeeze(ncv['agg_esd'][:,laysize])  #dsize #dsize is for dynamical
Gsize = squeeze(ncv['agg_G'][:,laysize])
G = squeeze(ncv['agg_G'][:,lay])
	#	ws = squeeze(ncv['agg_resuspensionlpm'][:,lay]) #(ncv['agg_ws'][:,lay])
		#ws2 = squeeze(ncv['agg_sinkinglpm'][:,lay]) #agg_Breakup
		#ws = -ws1 +ws2
aindex = np.diff(np.sign(np.diff(esd))).nonzero()[0] + 1             # local min & max
ymin = (np.diff(np.sign(np.diff(esd))) > 0).nonzero()[0] + 1         # local min
ymax = (np.diff(np.sign(np.diff(esd))) < 0).nonzero()[0] + 1         # local max# +1 
	#				axs[0].plot(secs/3600, totlpm,label=sys.argv[1]+ii+sys.argv[3]+jj) #+0.27 #SPMC	vs t
	#				axs[1].plot(secs/3600,esd*1e6,label=sys.argv[1]+ii+sys.argv[3]+jj) #+0.27#in m #ESD t
	#				axs[1].plot(secs[ymin]/3600, esd[ymin]*1e6, "o", label="min", color='r') #plot min
	#				axs[1].plot(secs[ymax]/3600, esd[ymax]*1e6, "o", label="max", color='b') #plot max
	#				axs[2].plot(Gsize,esd*1e6,alpha=0.3)  #totlpmsize/sqrt(Gsize) 	#ESD vs G
	#				axs[3].scatter(totlpm,esd*1e6, alpha=0.3) # 	G,totlpm,	#SPMC vs G
	#				axs[4].plot(secs/3600,ws,label=sys.argv[1]+ii+sys.argv[3]+jj)#mu=1.1e-3	#Ws vs t
ptime=secs[tslice]/3600.
array = ptime
					#get index for value/values to compare with spmc data/size data
for value in datatime: #s1datatimecorr: #[0:84]: #time in spmc data
	value = find_nearest (array, value)
			#print(find_nearest(array, value))
	Gidx.append(value) 	#append theoretical G index at certain t
			#Gidx=np.array(Gidx)
for values in s1datatime: 	#time in size data
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
			
for i1 in length:
	if i1<8:
		b=0
				#b=(dataspmc[i1]-totlpmi[i1])**2 #here #MSE absolute change
	else:
		b=sqrt(((dataspmc[i1]-totlpmi[i1])/dataspmc[i1])**2) #relative change
	sumb=sumb+b
			#test.append(b)
			#test01=sum(test)/len(test)
	mse=sumb/len(length)
			
for j1 in lengths:
			#d=(s1datasize[j1]-esdi[j1]*1e6)**2
	if j1<6:
		d=0
	else:
		d=sqrt(((s1datasize[j1]-esdi[j1]*1e6)/s1datasize[j1])**2)
	sumd=sumd+d
			#test1.append(d)
			#test11=sum(test1)/len(test1)
	mses=sumd/len(lengths)
			#secs[c]/3600, esd[c]*1e6
				
valley=(max(esd[ymax[1]],esd[ymax[2]])-esd[ymin[2]])/min(esd[ymax[1]],esd[ymax[2]])/0.76157 if len(ymin) >= 4 else 0 			#0.76157 is the valley ratio in data#*1e6
			#(max(esd[ymax])- max(esd[ymin]))*1e6 #.max#-esd[b].max
errorsum=mse+mses+abs(1-valley)
f=open('../relativeerror.txt', 'a+')
#		parts=re.split('(\d+)',i.name) 
#		f.write(parts[1]+" ") #" "
#		f.write(parts[3]+" ") #" "
#		f.write(parts[5]+" ") #" "
		#f.write(jj+" ")
		#f.write(kk+" ")
f.write(pathsub+" ") #"ID"+
f.write(str(mse)+" ")#"spmc"+
f.write(str(mses)+" ") #"size"
f.write(str(mse+mses)+" ") #"spmc+size"
f.write(str(valley)+" ") #"valley"+
f.write(str(len(ymin))+" ") #4 for valley true if run for 8 h "peaknr"
f.write(str(errorsum)+'\n')#'\n'"errortot"
f.flush()
f.close()
		#axs[0].scatter(datatime, totlpmi,label=sys.argv[1]+ii+sys.argv[3]+jj)   #scatter of modelled SPMC
		#axs[1].scatter(s1datatime,1e6*np.array(esdi),label=sys.argv[1]+ii+sys.argv[3]+jj)
		#print ("relative error SPMC=",mse, "relative error ESD=",mses, "sum=",mse+mses, "valley=",valley, "mins(4) = ",len(ymin),"errorsum=",errorsum) #"SumMSE=",mse+mses
	#print (filelist)
