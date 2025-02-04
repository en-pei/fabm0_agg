import matplotlib
from matplotlib import rc
import netCDF4
from pylab import *
import sys
import netcdftime
import numpy

rc('font',**{'family':'serif','serif':['FreeSans']})

def lognorm(x,mu,sig):
    return 1./(x*sig*sqrt(2*pi))*exp(-((log(x)-mu)**2)/(2*sig**2))

def logmedian(x,med,sig):
    return lognorm(x,log(med),sig)

def logmode(x,mode,sig):
    mu=log(mode)+sig**2
    return  1./(x*sig*sqrt(2*pi))*exp(-((log(x)-mu)**2)/(2*sig**2))

#data = asarray(numpy.loadtxt('paperplot/bale2002_ssc.dat',delimiter=', '))

ncfile=sys.argv[1]#'bale2002.nc'
nc=netCDF4.Dataset(ncfile)
ncv=nc.variables

utime=netcdftime.utime(ncv['time'].units)
time=utime.num2date(ncv['time'][:])
secs=ncv['time'][:]

lay=8
lpm = squeeze(ncv['spm_spm'][:,lay])
agglpm = squeeze(ncv['agg_agglpm'][:,lay])

por=0.8

when = int(float(sys.argv[2]))

por=0.2+0.6*agglpm[when]/(agglpm[when]+50.)
maxsize=50.+50.*agglpm[when]/(agglpm[when]+50.)
sig=0.75*agglpm[when]/(agglpm[when]+50.)

totlpm=lpm[when]+agglpm[when]/(1-por)


x=asarray([10**xx for xx in arange(0.0,3.5,0.01)])


plpm=logmode(x,20.0,0.8)
slpm=lpm[when]*plpm
pagg=logmode(x,maxsize,sig)
sagg=agglpm[when]/(1.0-por)*pagg

stot=slpm+sagg
pdf = stot/trapz(stot,x)
smean=trapz(pdf*x,x)

dpi=200
fig=figure(figsize=(5,5),dpi=dpi)
fig.subplots_adjust(left=0.15)

fill_between(x,pdf,0,color=(0.9,0.9,1.0))
semilogx(x,plpm*lpm[when]/(totlpm),'-',color='gray',lw=2.0)
semilogx(x,pagg*agglpm[when]/((1-por)*totlpm),'-',color='gray',lw=2.0)
semilogx(x,pdf,'k-',lw=2.0)
semilogx([smean,smean],[0,0.7*max(pdf)],'k--')
ylim(0,0.01)
xlim(1,10000)
text(smean,0.7*max(pdf),r'mean ESD = %0.1f $\mu m$'% smean)
xlabel(r'ESD [$\mu m$]')
ylabel('probability density')
title('size distribution at %0.1f h' % (secs[when]/3600.))

savefig('ESD_%s.png'%int(when),dpi=dpi)
