import matplotlib
from matplotlib import rc
from matplotlib import ticker, cm
import netCDF4
from pylab import *
import sys
import netcdftime
import numpy
#from astropy.visualization import simple_norm


def lognorm(x,mu,sig):
    return 1./(x*sig*sqrt(2*pi))*exp(-((log(x)-mu)**2)/(2*sig**2))

def logmedian(x,med,sig):
    return lognorm(x,log(med),sig)

def logmode(x,mode,sig):
    mu=log(mode)+sig**2
    return  1./(x*sig*sqrt(2*pi))*exp(-((log(x)-mu)**2)/(2*sig**2))

def size(lpm,agg):
    x=asarray([10**xx for xx in arange(0.5,3.5,0.01)])
    xnum=len(x)

    maxsize=50.+200*agglpm/(agglpm+50.)
    sig=0.75*agglpm/(agglpm+50.)

    tnum=len(lpm)
    pagg = zeros((tnum,xnum))
    slpm = zeros((tnum,xnum))
    sagg = slpm.copy()
    pdf = slpm.copy()
    smean = slpm.copy()

    for when in range(len(lpm)):
        plpm=logmode(x,50.0,0.5)
        pagg[when]=logmode(x,maxsize[when],sig[when])

        slpm[when]=lpm[when]*plpm
        sagg[when]=agglpm[when]*pagg[when]

        stot=slpm+sagg
        pdf[when] = stot[when]/trapz(stot[when],x)
        smean[when]=trapz(pdf[when]*x,x)
    return smean

rc('font',**{'family':'serif','serif':['FreeSans']})

#datadir = 'Plot/data'
data = asarray(numpy.loadtxt('SSC_20cmab_obs.dat',delimiter=', '))#'%s/SSC_20cmab_obs.dat'%datadir,
sizedata = asarray(numpy.loadtxt('D50_ESD_13cmab_obs.dat',delimiter=', ')) #'%s/D50_ESD_13cmab_obs.dat'%datadir,

ncfile='../../bale2002.nc'#'Results/v4h/results_vcurrent.nc' #sys.argv[1] #default output file
ncfile1='../../Results/v4hmoretep/results_vcurrent.nc'#'bale2002.nc'#
ncfile2='../../Results/v4hhalftep/results_vcurrent.nc' #'Results/vwithslackwater/results_vcurrent.nc'#'balemethod2.nc'

nc=netCDF4.Dataset(ncfile)
nc1=netCDF4.Dataset(ncfile1)
nc2=netCDF4.Dataset(ncfile2)

ncv=nc.variables
ncv1=nc1.variables
ncv2=nc2.variables

utime=netcdftime.utime(ncv['time'].units)
time=utime.num2date(ncv['time'][:])
secs=ncv['time'][:]

lay=16 #3 #16 # 11 # layer for SPM
lpm = squeeze(ncv['spm_spm'][:,lay])
agglpm = squeeze(ncv['agg_agglpm'][:,lay])
totlpm=lpm+agglpm

lpm1 = squeeze(ncv1['spm_spm'][:,lay])
agglpm1 = squeeze(ncv1['agg_agglpm'][:,lay])
totlpm1=lpm1+agglpm1

lpm2 = squeeze(ncv2['spm_spm'][:,lay])
agglpm2 = squeeze(ncv2['agg_agglpm'][:,lay])
totlpm2=lpm2+agglpm2

h=3600./60 #/120. output was per 120s
tstart=int(0*h)#int(6.*h)
tend=tstart+int(39*h)#tstart+int(8*h)
tslice=slice(tstart,tend)
ptime=secs[tslice]/3600.

laysize=9
Gsize = squeeze(ncv['agg_G'][:,laysize])
G = squeeze(ncv['agg_G'][:,lay])
ws = squeeze(ncv['agg_ws'][:,lay])

lpm1 = squeeze(ncv1['spm_spm'][:,lay])
agglpm1 = squeeze(ncv1['agg_agglpm'][:,lay])
G1 = squeeze(ncv1['agg_G'][:,lay])
#ws1 = squeeze(ncv1['agg_ws'][:,lay])

G2 = squeeze(ncv2['agg_G'][:,lay])

f=figure(figsize=(9,2.5))  #plotting G
#f.subplots_adjust(bottom=0.15)
lw=2.0
plot(ptime,G[tslice],'k-',lw=lw, label='noslackwater')#'2x TEP'
plot(ptime,Gsize[tslice],'r--',lw=lw, label='G0.13')
#plot(ptime,G2[tslice],'k--',lw=lw, label='withslackwater')
xlabel('time (h)')
ylabel(r'turbulence shear (1/s)')
plt.xlim(0,)
#legend(loc='upper right') 
savefig('G.png'); close()


esd=squeeze(ncv['agg_esd'][:,laysize])

s1datatime=sizedata[:,0]-5
s1datasize=sizedata[:,1]


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx#, array[idx]

Gidx=[]
array = ptime
for value in s1datatime[0:84]:
	value = find_nearest (array, value)
	#print(find_nearest(array, value))
	Gidx.append(value)
Gidx=np.array(Gidx)
#print Gidx
Gmodel=[]
Dmodel=[]

for i in Gidx:
	Gi=G[i]
	Di=esd[i+17] #+0.28h*60
	#Di=esd[i]
	Gmodel.append(Gi)
	Dmodel.append(Di)
#print Gmodel, Dmodel

Gmodel=np.array(Gmodel)	
Dmodel=np.array(Dmodel)	

fig,ax=plt.subplots(1,1) #plot ESD from model and data against G shear
lw=2.0
ax.scatter(Gmodel[5:],Dmodel[5:]*1e6*0.6+50,label='simulated') #Dmodel[1:]*1e6*0.6+50,
ax.scatter(Gmodel[5:],s1datasize[5:84],label='data',marker='^') #s1datasize[0:84],
ax.set_xlabel('turbulence shear (1/s)')
#ax.set_yscale('log')
ax.set_ylabel('ESD ($\mu$m)')
ax.legend()
#savefig('ESDvsG.png'); close()

print s1datatime, s1datasize, Gmodel, Dmodel

f=figure(figsize=(9,2.5)) 
plot(ptime[5:]+0.28,esd[5:]*1e6*0.6+50,label='simulated') #ptime[1:]+0.28,esd[1:]*1e6*0.3+70,
plt.scatter(s1datatime[5:],s1datasize[5:],label='data')
xlabel('time (h)')
ylabel('ESD ($\mu$m)')
plt.xlim(0,)
plt.legend()
plt.show()

#legend(loc='upper right') 
#savefig('G.png'); close()
