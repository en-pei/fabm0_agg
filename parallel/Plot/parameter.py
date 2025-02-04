import matplotlib
from matplotlib import rc
from matplotlib import ticker, cm
import netCDF4
from pylab import *
from scipy.interpolate import interp1d
import netcdftime
import numpy	
import os
import sys
import subprocess
#from glob import glob

#get data from bale2002 observation
datadir = '/home/enpei/mossco/setups/bale2002/Plot/data'
data = asarray(numpy.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=', '))
sizedata = asarray(numpy.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))
s1datatime=sizedata[:,0]-5
s1datasize=sizedata[:,1] #data is in micrometer
datatime=data[:,0]-5.
dataspmc=data[:,2]

lay= 16 #layer for SPM concentration
laysize= 9 #layer for ESD



#glob("/home/enpei/mossco/setups/bale2002/paramstudy/*/")
#subfolders = [ f.name for f in os.scandir('/home/enpei/mossco/setups/bale2002/paramstudy') if f.is_dir() ]
#print (subfolders)


filelist=sys.argv[2].split() #split the string
filelist1=sys.argv[4].split()
filelist2=sys.argv[5].split()
#print sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]
#currentDirectory = os.getcwd()
#print currentDirectory
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx#, array[idx]

h=3600./60 #/120. output was per 120s
tstart=int(0*h)#int(6.*h)
tend=tstart+int(39*h)#tstart+int(8*h)
tslice=slice(tstart,tend)

Gidx=[]
Gsidx=[]



s1datatimecorr=list(i for i in s1datatime) #(i-0.27) #time in datasize, 2.7 h earlier
#print len(datatime) #79
length=range(len(datatime)) 				#idexing spmc data
lengths=range(len(s1datasize)) 				#indexing size data
#print length, lengths


fig, axs = plt.subplots(5, sharex=False)
for ii in filelist:
#	if os.path.isdir(sys.argv[1]+ii): 
#		print (sys.argv[1],ii,'existing')
#		continue	
#	except IOError: #(IOError, OSError): #except FileNotFoundError:
#		continue
	for jj in filelist1:
		for kk in filelist2:
			#print (sys.argv[1],'are',ii,sys.argv[3],'are',jj,sys.argv[6],'are',kk)
			#if os.path.isdir(sys.argv[1]+ii+sys.argv[3]+jj+sys.argv[6]+kk):
			#	continue
			if os.path.isdir(sys.argv[1]+ii+sys.argv[3]+jj+sys.argv[6]+kk):
				#print ("sucess")
			#try:
				ncfile=sys.argv[1]+ii+sys.argv[3]+jj+sys.argv[6]+kk+'/bale2002.nc' #find nc 
			#except IOError: #(IOError, OSError): #except FileNotFoundError:
			#	print (sys.argv[1],ii,sys.argv[3],jj,sys.argv[6],kk, 'No nc file')
			#	continue
#		else:
#			continue
				print 'ncfiles are',ncfile
				nc=netCDF4.Dataset(ncfile)
				ncv=nc.variables
				utime=netcdftime.utime(ncv['time'].units)
				time=utime.num2date(ncv['time'][:])
				secs=ncv['time'][:]
				#print len(secs)
				if len(secs)<481:
					continue #skip run time smaller ones
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
				ws = squeeze(ncv['agg_resuspensionlpm'][:,lay]) #(ncv['agg_ws'][:,lay])
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
				#	test.append(b)
				#test01=sum(test)/len(test)
				mse=sumb/len(length)
			
				for j1 in lengths:
				#d=(s1datasize[j1]-esdi[j1]*1e6)**2
					if j1<6:
						d=0
					else:
						d=sqrt(((s1datasize[j1]-esdi[j1]*1e6)/s1datasize[j1])**2)
					sumd=sumd+d
				#	test1.append(d)
				#test11=sum(test1)/len(test1)
				mses=sumd/len(lengths)
				#secs[c]/3600, esd[c]*1e6
			
				valley=(max(esd[ymax[1]],esd[ymax[2]])-esd[ymin[2]])/min(esd[ymax[1]],esd[ymax[2]])/0.76157 if len(ymin) >= 4 else 0 			#0.76157 is the valley ratio in data#*1e6
				#(max(esd[ymax])- max(esd[ymin]))*1e6 #.max#-esd[b].max
				errorsum=mse+mses+(1-valley)
				f=open('relativeerror.txt', 'a+')
				f.write(ii+" ") #" "
				f.write(jj+" ")
				f.write(kk+" ")
				f.write(str(mse)+" ")
				f.write(str(mses)+" ")
				f.write(str(mse+mses)+" ")
				f.write(str(valley)+" ")
				f.write(str(len(ymin))+" ") #4 for valley true if run for 8 h
				f.write(str(errorsum)+'\n')#'\n'
			#axs[0].scatter(datatime, totlpmi,label=sys.argv[1]+ii+sys.argv[3]+jj)   #scatter of modelled SPMC
			#axs[1].scatter(s1datatime,1e6*np.array(esdi),label=sys.argv[1]+ii+sys.argv[3]+jj)
				print "relative error SPMC=",mse, "relative error ESD=",mses, "sum=",mse+mses, "valley=",valley, "mins(4) = ",len(ymin),"errorsum=",errorsum #"SumMSE=",mse+mses
				os.remove(ncfile)
				os.remove(os.path.isdir(sys.argv[1]+ii+sys.argv[3]+jj+sys.argv[6]+kk/*)
				#axs[4].scatter(dataspmc[8:],totlpmi[8:]) #plot modelled with observation SPMC
				esdi=np.array(esdi)
			#axs[4].scatter(s1datasize[6:],esdi[6:]*1e6)  #data vs model
		#	msei.append(mse)
		#	msesi.append(mses)
		#nmse=msei/(sum(dataspmc)/len(length)) 						#normalized mean squre error
		#nmsei=msesi/(sum(s1datasize)/len(lengths))
		#print "NMSE of SPMC=",msei, "NMSE of ESD=",msesi,"MSE of SPMC=",nmse,"MSE of ESD=",nmsei
			else: 
				print (sys.argv[1],ii,sys.argv[3],jj,sys.argv[6],kk,"fail"))				
				#pass
				#os.remove(ncfile)
				continue
		#else: 
		#	continue
axs[0].plot(datatime[8:],dataspmc[8:],'k^',ms=5,label='total lit, data') 		#dataSPMCvs	t
axs[0].set_ylabel('SPMC (g m$^{-3}$)')
#axs[0].set_xlim(0,16)
axs[1].plot(s1datatime[6:],s1datasize[6:],'g^',ms=5.,label='ESD, data') 			#dataESD vs	t

#dataindex = np.diff(np.sign(np.diff(s1datasize[6:]))).nonzero()[0] + 1+6               # local min & max
#datamin = (np.diff(np.sign(np.diff(s1datasize[6:]))) > 0).nonzero()[0] + 1+6         # local min
#datamax = (np.diff(np.sign(np.diff(s1datasize[6:]))) < 0).nonzero()[0] + 1+6         # local max# +1 due to
#datavalley=(min(s1datasize[datamax])- min(s1datasize[datamin]))
#print "datavalley=",datavalley


axs[1].set_xlabel('time (h)')
axs[1].set_ylabel('mean ESD ($\mu m$)')
#axs[1].set_xlim(0,16)
#axs[4].plot(s1datasize[6:],s1datasize[6:],'k')
Gmodel=[]
Dmodel=[]
for i in Gsidx:
	Gi=Gsize[i] #Gsize at certain time for plotting D vs G
	#Di=esd[i+17] #esd[i+17]+0.28h*60
	#Di=esd[i]
	#Dmodel.append(Di)
	Gmodel.append(Gi)
#Gmodel=np.array(Gmodel)	
#Dmodel=np.array(Dmodel)
#print len(Gmodel) #87 ?261
#print len(s1datasize)	#87
#axs[2].scatter(Gmodel,s1datasize,label='data') #[5:84]  				#dataESD vs	G
axs[2].set_xlabel('G')
axs[2].set_ylabel('ESD ($\mu$m)')
axs[3].set_xlabel('SPMC (g m$^{-3}$)')
axs[3].set_ylabel('ESD ($\mu$m)')
axs[4].set_ylabel('dD (m s$^{-1}$)') #'ws (m s$^{-1}$)' ('modelled SPMC')# ('modelled ESD ($\mu$m)')
#plt.legend()
#plt.show()



###################################################################################################################################

#fig,ax=plt.subplots(1,1) #plot ESD from model and data against G shear
#lw=2.0
#ax.scatter(Gmodel[5:],Dmodel[5:]*1e6*0.6+50,label='simulated') #Dmodel[1:]*1e6*0.6+50,
#ax.scatter(Gmodel[5:],s1datasize[5:84],label='data',marker='^') #s1datasize[0:84],
#ax.set_xlabel('turbulence shear (1/s)')
##ax.set_yscale('log')
#ax.set_ylabel('ESD ($\mu$m)')
#ax.legend()
#savefig('ESDvsG.png'); close()

#print s1datatime, s1datasize, Gmodel, Dmodel

#f=figure(figsize=(9,2.5)) 
#plot(ptime[5:]+0.28,esd[5:]*1e6*0.6+50,label='simulated') #ptime[1:]+0.28,esd[1:]*1e6*0.3+70,
#plt.scatter(s1datatime[5:],s1datasize[5:],label='data')
#xlabel('time (h)')
#ylabel('ESD ($\mu$m)')
#plt.xlim(0,)
#axs.legend()
#plt.show()

#CONTAINED IN THE LOOP
#ncfile1='paramstudy/test.nc'#'Results/v4hmoretep/results_vcurrent.nc'#'bale2002.nc'#
#ncfile2='Results/v4hhalftep/results_vcurrent.nc' #'Results/vwithslackwater/results_vcurrent.nc'#'balemethod2.nc'

#nc=netCDF4.Dataset(ncfile)
#nc1=netCDF4.Dataset(ncfile1)
#nc2=netCDF4.Dataset(ncfile2)

#ncv=nc.variables
#ncv1=nc1.variables
#ncv2=nc2.variables


#lay=16 #layer for SPM concentration
#lpm = squeeze(ncv['spm_spm'][:,lay])
#agglpm = squeeze(ncv['agg_agglpm'][:,lay])
#totlpm=lpm+agglpm

#lpm1 = squeeze(ncv1['spm_spm'][:,lay])
#agglpm1 = squeeze(ncv1['agg_agglpm'][:,lay])
#totlpm1=lpm1+agglpm1



##interpolation for calculating the standard deviation 
#datalpmtime=data[:,0]-5
#xarray=np.asarray(datalpmtime)
#def find_nearest(array,num):
#	nearest_idx = np.where(abs(xarray-num)==abs(xarray-num).min())[0]
#	nearest_val = array[abs(xarray-num)==abs(xarray-num).min()] 
#	idxmin = nearest_idx-1
#	idxmax = nearest_idx+1
#	valmin = xarray[idxmin]
#	valmax = xarray[idxmax]
#	variablemin = array[idxmin]
#	variablemax = array[idxmax]
#	valx=np.array([idxmin,idxmax]).squeeze()
#	valy=np.array([valmin,valmax]).squeeze()
#	y_interp = interp1d(valy, valx)
#	return (num, y_inter(num))

#find_nearest(totlpm,10) #here!

#print valx.shape, valy.shape, valx, valy
#print nearest_idx, nearest_val, idxmin, valmin, idxmax, valmax


