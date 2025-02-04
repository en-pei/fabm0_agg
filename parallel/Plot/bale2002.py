import matplotlib
from matplotlib import rc
from matplotlib import ticker, cm
import netCDF4
from pylab import *
import sys
import netcdftime
import numpy
import pickle
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



datadir = 'Plot/data'
data = asarray(numpy.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=', '))
sizedata = asarray(numpy.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))

with open ('myoutput.txt', 'rb') as fp:
	itemlist = pickle.load(fp)
#print(itemlist[0])

ncfile='bale2002.nc'#'Results/v4h/results_vcurrent.nc' #sys.argv[1] #default output file
ncfile1='Results/v4hmoretep/results_vcurrent.nc'#'bale2002.nc'#
ncfile2='Results/v4hhalftep/results_vcurrent.nc' #'Results/vwithslackwater/results_vcurrent.nc'#'balemethod2.nc'

nc=netCDF4.Dataset(ncfile)
nc1=netCDF4.Dataset(ncfile1)
nc2=netCDF4.Dataset(ncfile2)

ncv=nc.variables
ncv1=nc1.variables
ncv2=nc2.variables

utime=netcdftime.utime(ncv['time'].units)
time=utime.num2date(ncv['time'][:])
secs=ncv['time'][:]

lay=16 #3 #16 # 11 # layer for SPMC
laysize=9 #10 #6 # 7 # layer for size


lpm = squeeze(ncv['spm_spm'][:,lay])
agglpm = squeeze(ncv['agg_agglpm'][:,lay])
totlpm=lpm+agglpm

#size layer SPMC
lpm1 = squeeze(ncv['spm_spm'][:,laysize])
agglpm1 = squeeze(ncv['agg_agglpm'][:,laysize]) #attention to the ncv
totlpm1=lpm1+agglpm1

lpm2 = squeeze(ncv2['spm_spm'][:,lay])
agglpm2 = squeeze(ncv2['agg_agglpm'][:,lay])
totlpm2=lpm2+agglpm2

h=3600./60 #/120. output was per 120s
tstart=int(6*h)#int(6.*h) #inital time for plotting
tend=tstart+int(39*h)#tstart+int(8*h)
tslice=slice(tstart,tend)
ptime=secs[tslice]/3600.


dpi=200
figure(figsize=(9,9*0.75))
lw=2.0
plot(ptime,totlpm[tslice],'k-',lw=lw+1,label='simulated SPM' )#1xTEP'total lit3''TSM')
plot(ptime,totlpm1[tslice],'r--',lw=lw+1,label= 'size layer SPMC')#'total lit2''TSM') #vcurrent
#plot(ptime,totlpm2[tslice]*0.001,'b:',lw=lw+1,label= '0.5xTEP')#'total lit2''TSM')

#plot(ptime,agglpm[tslice],'r-',lw=lw,label='aggregates mass')

#plot(ptime,agglpm1[tslice]*0.001,'r--',lw=lw,label='lit in aggregates1')
#plot(ptime,agglpm2[tslice]*0.001,'r--',lw=lw,label='lit in aggregates2')
#plot(ptime,agglpm[tslice]/totlpm[tslice]*100,label='aggregates in total mass (%)')

#plot(ptime,lpm[tslice]*0.001,'g-',lw=lw,label='fine SPM')#'fine SPM')

#plot(ptime,lpm2[tslice]*0.001,'g-',lw=lw,label='lit in water2')#'fine SPM')

plot(data[:,0]-5,data[:,2],'k^',ms=5.,label='SPM, data')#'TSM\nBale.ea2002')
#plot(ptime,totlpm[tslice]/totlpm[tslice]+3500, color='r', linestyle='-') #add a line to see if mass conserves
tick_params(axis='both', which='minor', labelsize=20)
legend(loc='upper right') 
#legend(loc='center left', bbox_to_anchor=(1, 0.5))
xlabel('time (h)', fontsize=18)
ylabel('lithogenous concentration (g m$^{-1}$)', fontsize=18)
#ylim(0,4.5)
xlim(6,14)
#xlim(0,40)
savefig('SSC.png',dpi=dpi); close()
totlpm=np.array(totlpm)
totlpm1=np.array(totlpm1)
#print (min(totlpm[200:]),min(totlpm1[200:]))



##vertical profile, SPMC g/m-3, log scale
dpi=400
##construct array for contour plot
agglpmz = squeeze(ncv['agg_agglpm'][:,:])
agglpmz = transpose(agglpmz)
lpmz = squeeze(ncv['spm_spm'][:,:])
z=squeeze(ncv['z'][:,:])
aggmassz=squeeze(ncv['agg_aggmass'][:,:])
aggmassz = transpose(aggmassz)
lpmz = transpose(lpmz)
spmcz=aggmassz+lpmz
timez=np.c_[secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs]
#timez=np.c_[secs,secs,secs,secs,secs,secs]
timez=transpose(timez)
datatime=np.c_[data[:,0]-5,data[:,0]-5,data[:,0]-5]
dataz=np.c_[-(0.28-data[:,1]*0.01),-(0.28-data[:,1]*0.01),-(0.28-data[:,1]*0.01)] #from cm to m #-(0.28-23*0.01)
dataspmc=np.c_[data[:,2],data[:,2],data[:,2]]
fig,ax=plt.subplots(1,1)
#cp = ax.contourf(timez/3600,transpose(z),aggmassz,locator=ticker.LogLocator(),cmap=cm.OrRd)
#cp = ax.contourf(timez/3600,transpose(z),lpmz,locator=ticker.LogLocator(),cmap=cm.OrRd)
#cp = ax.contourf(timez/3600,transpose(z),spmcz,locator=ticker.LogLocator(),cmap=cm.RdBu_r) #norm=matplotlib.colors.SymLogNorm(linthresh=1e6, linscale=3000, vmin=spmcz.min(),
		#vmax=spmcz.max()), vmin=spmcz.min(), vmax=spmcz.max(),
#cbar=fig.colorbar(cp)# # Add a colorbar to a plot
cp = ax.contourf(timez/3600,transpose(z),np.log10(spmcz),cmap='Reds') #np.log10(spmcz) #norm=matplotlib.colors.LogNorm()(vmin=spmcz.min(), vmax=spmcz.max()))
#norm=matplotlib.colors.PowerNorm(gamma=1. / 2.),
#                       cmap='PuBu_r')
colorbar=fig.colorbar(cp, ax=ax, extend='max')
colorbar.set_label('SPMC',labelpad=-40, y=1.05, rotation=0) #log(SPMC)
##other stuff
#ctks = [0,1e2,10e2,50e2,100e2,500e2,1000e2,5000e2]
#cbar.set_ticks(ctks)
#print (lpmz[0,:])

ax.set_xticks(ax.get_xticks()[::1])
ax.set_xlabel('time (h)')
ax.set_ylabel('water depth (m)')
#im=plt.scatter(datatime,dataz,c=dataspmc,cmap=cm.OrRd,marker='^',s=5.) #c=log10(dataspmc), vmin=aggmassz.min(), vmax=aggmassz.max())
#plt.clim(-3,6)
#plt.clim(2000,3600)
#plt.clim(0,56000)
#fig.colorbar(im)
savefig('verticalSSC.png'); close()
#plt.show()


##vertical G
zesd = squeeze(ncv['z'][:,:])
zesd = transpose(zesd)
gd = squeeze(ncv['agg_G'][:,:]) #testing u component#(ncv['u'][:,:])
gd = transpose(gd)
fig,ax=plt.subplots(1,1)
ax.set_xlabel('time (h)')
ax.set_ylabel('water depth (m)')
cpg = ax.contourf(timez/3600,zesd,gd,cmap=cm.Blues)#log(gd)#,locator=ticker.LogLocator())
#levels=[0,0.001,0.01,0.1,1,5,10]
#norm = simple_norm(image, 'sqrt')
clbg=fig.colorbar(cpg, ax=ax)#,origin='lower',norm=norm)#,ticks=levels)
clbg.set_label('G',labelpad=-40, y=1.05, rotation=0) #u component
savefig('gvertical.png',dpi=dpi); close
#plt.show()

#G profile at certain time
fig,ax11=plt.subplots(1,1)
gd1 = transpose(gd)
zesd1 = transpose(zesd)
#ax11.plot(gd1[1],zesd1[0],label='0 h')
#ax11.plot(gd1[60],zesd1[0],label='1 h')
ax11.plot(gd1[60*2],zesd1[0],label='2 h')
ax11.plot(gd1[60*3],zesd1[0],label='3 h')
ax11.plot(gd1[60*4],zesd1[0],label='4 h')
#ax11.plot(gd1[60*4],zesd1[0],label='5 h')
ax11.set_xlabel('G (s$^{-1}$)')
ax11.set_ylabel('water depth (m)')
ax11.legend(loc='upper right')
savefig('gprofile.png',dpi=dpi)#; close
#print (len(gd),len(zesd1))
print (np.max(gd))

##vertical Dsize
Dsize = squeeze(ncv['agg_Dsize'][:,:])
Dsize = transpose(Dsize)
fig,ax=plt.subplots(1,1)
ax.set_xlabel('time (h)')
ax.set_ylabel('water depth (m)')
cpg = ax.contourf(timez/3600,zesd,Dsize*1e6,cmap=cm.coolwarm)#,locator=ticker.LogLocator())
#levels=[0,0.001,0.01,0.1,1,5,10]
#norm = simple_norm(image, 'sqrt')
clbg=fig.colorbar(cpg, ax=ax)#,origin='lower',norm=norm)#,ticks=levels)
clbg.set_label('Dsize($\mu$m)',labelpad=-40, y=1.05, rotation=0)
savefig('Dsizevertical.png',dpi=dpi); close



fig,ax2=plt.subplots(1,1)
cp2 = ax2.contourf(timez/3600,transpose(z),aggmassz/(spmcz),cmap=cm.OrRd) #locator=ticker.LogLocator(),
ax2.set_xticks(ax.get_xticks()[::1])
ax2.set_xlabel('time (h)')
ax2.set_ylabel('water depth (m)')
cbar2=fig.colorbar(cp2)
#plt.show()
savefig('verticalaggratio.png'); close()

G = squeeze(ncv['agg_G'][:,lay])
ws = squeeze(ncv['agg_ws'][:,lay])

#changed to size layer
lpm1 = squeeze(ncv['spm_spm'][:,laysize])
agglpm1 = squeeze(ncv['agg_agglpm'][:,laysize])
G1 = squeeze(ncv['agg_G'][:,laysize])
ws1 = squeeze(ncv['agg_ws'][:,laysize])

Gesd = squeeze(ncv['agg_G'][:,laysize]) #size layer turbulence

f=figure(figsize=(9,2.5),dpi=dpi)  #plotting G
#f.subplots_adjust(bottom=0.15)
lw=2.0
plot(ptime,G[tslice],'k-',lw=lw, label='G SPMC layer')#'2x TEP'
plot(ptime,Gesd[tslice],'k--',lw=lw, label='G size layer')
xlabel('time (h)')
ylabel(r'turbulence shear (1/s)')
plt.xlim(0,)
legend(loc='upper right') 
savefig('G.png',dpi=dpi); close()
G=np.array(G)
Gesd=np.array(Gesd)

print (max(Gesd[200:]), max(G[200:]))

## 2 axis plot G & SSC
#fig, ax1 = plt.subplots()
##color = 'tab:k'
#ax1.set_xlabel('time [h]')
#ax1.set_ylabel('SSC (g/L)')#, color=color
#ax1.plot(ptime, totlpm[tslice],color='k')
#ax1.plot(data[:,0]-5.0,data[:,1],'k^',ms=5.,label='total lit, data')#'TSM\nBale.ea2002')
#ax1.tick_params(axis='SSC (g/m$^{-3}$)')#, labelcolor=color

#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#color = 'tab:red'
#ax2.set_ylabel('turbulence shear (1/s)', color=color)  # we already handled the x-label with ax1
#ax2.plot(ptime, G[tslice],'--',color=color)
#ax2.tick_params(axis='turbulence shear (1/s)', labelcolor=color)
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#savefig('G.png',dpi=dpi); close()


lpm = squeeze(ncv['spm_spm'][:,laysize])
agglpm = squeeze(ncv['agg_agglpm'][:,laysize])
agglpm1 = squeeze(ncv1['agg_agglpm'][:,laysize])
agglpm2 = squeeze(ncv2['agg_agglpm'][:,laysize])
G = squeeze(ncv['agg_G'][:,laysize])
G1 = squeeze(ncv['agg_G'][:,laysize])
#G2 = squeeze(ncv2['agg_G'][:,laysize])

figure(figsize=(9,4.5),dpi=dpi)
#f.subplots_adjust(bottom=0.15)

#lw=2.0
##plot(ptime,size(lpm[tslice],agglpm[tslice]),'k-',lw=lw)
#plot(ptime,0.01+0.02*agglpm[tslice]/G[tslice],'k-',lw=lw,label='2x TEP')
##plot(ptime,0.01+0.02*agglpm1[tslice]/G1[tslice],'k--',lw=lw,label='1x TEP') #lpm1 is with kc=11, while the other one is kc=22
#xlabel('time [h]')
#ylabel(r'mean ESD [$\mu m$]')
#legend()
#savefig('Size.png',dpi=dpi); close()


#2 axis plot ESD & G
esd=squeeze(ncv['agg_esd'][:,laysize])
esd1=squeeze(ncv1['agg_esd'][:,laysize])
esd2=squeeze(ncv2['agg_esd'][:,laysize])
dsizel=squeeze(ncv['agg_Dsize'][:,laysize]) #Dsize for sizelayer
##subplots 2 ssc and size
fig, (ax1, ax3) = plt.subplots(2, sharex=True)
#fig, ax1 = plt.subplots()
#color = 'tab:k'
#ax1.set_xlabel('time (h)')
ax1.set_ylabel('size $\mu$m')#, color=color

#ax1.plot(ptime, 0.01+0.02*agglpm[tslice]/G[tslice],color='k') #old formulation
#ax1.plot(ptime, 0.01+0.02*agglpm2[tslice]/G2[tslice],'--',color='k') #old formulation

##method3 #0.0001_rk + 0.0003_rk*1.d-3*aggmass/G  
#ax1.plot(ptime, 0.1*1e6*(0.0001+0.0003*1e-3*agglpm[tslice]/G[tslice]),color='k') #method 3 formulation
#ax1.plot(ptime, 0.1*1e6*(0.0001+0.0003*1e-3*agglpm2[tslice]/G2[tslice]),'--',color='k') #method 3 formulation

##directly using model output

##here, 0D and 1D diagnostic 
#ax1.plot(ptime, 1e6*esd[tslice],color='green',lw=lw+1,label='1D diagnostic $D$',linestyle='--') #1xTEPesd model output#
#ax1.plot(itemlist[0],itemlist[1],label='0D simulated $D$',color='green',lw=lw+1)
ax1.plot(ptime,1e6*dsizel[tslice],label='1D dynamical $D$',color='red',lw=lw+1)

#ax1.plot(ptime, 1e6*esd1[tslice],'--',color='red',label='1.5xTEP') #esd model output
#ax1.plot(ptime, 1e6*esd2[tslice],':',color='blue',label='0.5xTEP') #esd model output
#ax1.set_ylim(0,300)
ax1.legend(loc='upper right')

ax3.plot(ptime,totlpm[tslice],'k--',lw=lw+1,label='1D diagnostic $D$')#1xTEP'total lit3''TSM')
ax3.plot(data[:,0]-5,data[:,2],'k^',ms=5.,label='data')#'TSM\nBale.ea2002')
#ax3.plot(itemlist[0],[item*1000 for item in itemlist[2]],label='0D')  #wrong depth shouldn't be here...
#tick_params(axis='both', which='minor', labelsize=20)
#legend(loc='upper right') 
#legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax3.set_xlabel('time (h)')
ax3.set_ylabel('SPM concentration (g m$^{-1}$)')
ylim(0,4500)
#xlim(6,14)
#xlim(0,40)

s1datatime=sizedata[:,0]-5
s1datasize=sizedata[:,1] #data is in micrometer
ax1.scatter(s1datatime, s1datasize,color='green',label='datasize')

ax1.set_ylim(0,)
#ax1.tick_params(axis='size')#, labelcolor=color

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax2.set_ylabel('turbulence shear (1/s)', color=color)  # we already handled the x-label with ax1
ax2.plot(ptime, G[tslice],linewidth=3, color=color,alpha=0.2)
#ax2.tick_params(axis='turbulence shear (1/s)', labelcolor=color)
#ax1.set_ylim(0,400)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.show()
savefig('Size.png',dpi=dpi); close()



fig,ax=plt.subplots(1,1)
ax.set_xlabel('G(1/s)') #time (h)
ax.set_ylabel('size change term (m/s)')
coagulationlpm=squeeze(ncv['agg_coagulationlpm'][:,laysize]) 
coagulationcalculated=14*0.5*G1[tslice]*lpm1[tslice]*1e-3*(dsizel[tslice])**2
diffsetlpm=squeeze(ncv['agg_diffsetlpm'][:,laysize]) 
sinkinglpm=squeeze(ncv['agg_sinkinglpm'][:,laysize]) 
###coagulation*1.d-3*lpm*Dsize**(4-self%fractal_dimension)
breakuplpm=squeeze(ncv['agg_Breakup'][:,laysize])
breakupcalculated=13000*(G1[tslice])**1.5*(dsizel[tslice]-1e-5)*(dsizel[tslice])**2
###-self%breakup_factor*G**1.5d0*(Dsize-self%min_size)*Dsize**2 
#ax.plot(ptime,1e-3*agglpm[tslice]*breakuplpm[tslice],'--',lw=lw+1,label='breakup rate') #breakup rate in kg/m-3/s
#ax.plot(ptime,1e-3*coagulationlpm[tslice],'-,',lw=lw+1,label='coagulation rate') #ptime G  #coagulation rate in kg/m-3/s
#ax.plot(ptime,-breakupcalculated,'--',lw=lw+1,label='breakup rate') #breakuplpm[tslice] , #breakup rate dynamical size change  #ptime  #G1 is laysize G
#ax.plot(ptime,coagulationcalculated,'-,',lw=lw+1,label='coagulation rate') #coagulationlpm[tslice] coagulation rate dynamical size change

ax.plot(ptime,-breakuplpm[tslice],'--',lw=lw+1,label='breakup rate') #G1[tslice]**1.5
ax.plot(ptime,coagulationlpm[tslice],'-,',lw=lw+1,label='coagulation rate') #totlpm1[tslice]
ax.plot(ptime,diffsetlpm[tslice],'-,',lw=lw+1,label='differential settling') #totlpm1[tslice]
ax.plot(ptime,sinkinglpm[tslice],'-,',lw=lw+1,label='sinking') #totlpm1[tslice]

legend()
savefig('coagulationbreakup.png',dpi=dpi); close
#print (breakuplpm[tslice]/coagulationlpm[tslice]) #*(G1[tslice])**0.5
#print (breakupcalculated/coagulationcalculated)


fig,ax=plt.subplots(1,1)
ax.set_xlabel('G(1/s)') #time (h)
ax.set_ylabel('size change term (m/s)')
dDsum=coagulationlpm[tslice]-breakuplpm[tslice]+diffsetlpm[tslice]+sinkinglpm[tslice]
ax.plot(G1[tslice],dDsum,label='total dD')
legend()
savefig('totaldD.png',dpi=dpi); close



f=figure(figsize=(9,4.5),dpi=dpi)
#f.subplots_adjust(bottom=0.15)
lw=2.0
#plot(ptime,size(lpm[tslice],agglpm[tslice]),'k-',lw=lw)
plot(ptime,1000*ws[tslice],'k-',lw=lw, label='size layer sinking')
plot(ptime,1000*ws1[tslice],'k--',lw=lw, label='SPMC layer sinking')
xlabel('time [h]')
ylabel(r'sinking velocity [mm/s]')
#ylim(-20,0)
legend()
savefig('Ws.png',dpi=dpi); close()
ws=np.array(ws)
ws1=np.array(ws1)
#print (max(ws[500:]),max(ws1[500:]))


##vertical profile, sinking velocity
dpi=200
##construct array for contour plot
wsz = squeeze(ncv['agg_ws'][:,:])
wsz = transpose(-wsz)  #????? wy the minus
fig,ax=plt.subplots(1,1)
cpw = ax.contourf(timez/3600,transpose(z),wsz,cmap=cm.YlOrBr) #locator=ticker.LogLocator(),
clws=fig.colorbar(cpw, ax=ax)
clws.set_label('$w_s$ (m s$^{-1}$)',labelpad=-40, y=1.05, rotation=0)
#plt.clim(0,1050)
#ax.set_xticks(ax.get_xticks()[::1])
#plt.show()
ax.set_xlabel('time (h)')
ax.set_ylabel('water depth (m)')
savefig('Wsvertical.png',dpi=dpi); close



##vertical profile, size, log scale
##D50_ESD_13cmab_obs.dat
dpi=200
##construct array for contour plot
sizez = squeeze(ncv['agg_esd'][:,:])
sizez = transpose(sizez*1e6)
sdatatime=np.c_[sizedata[:,0]-5,sizedata[:,0]-5,sizedata[:,0]-5]
sdataz=np.c_[-(0.28-sizedata[:,2]),-(0.28-sizedata[:,2]),-(0.28-sizedata[:,2])]  #-(0.28-0.13)
sdatasize=np.c_[sizedata[:,1],sizedata[:,1],sizedata[:,1]]

fig,ax=plt.subplots(1,1)
cps = ax.contourf(timez/3600,zesd,sizez,cmap=cm.Greens) #np.log10(sizez) locator=ticker.LogLocator(),
#ims=plt.scatter(sdatatime,sdataz,c=sdatasize,cmap=cm.PuBu,marker='^',s=5.)#, vmin=1e-6, vmax=1e-1)
clbsize=fig.colorbar(cps, ax=ax)
clbsize.set_label('($D$)',labelpad=-40, y=1.05, rotation=0) #log
#fig.colorbar(ims, ax=ax)
#plt.clim(0,1050)
ax.set_xticks(ax.get_xticks()[::1])
#plt.show()
ax.set_xlabel('time (h)')
ax.set_ylabel('water depth (m)')
savefig('Sizevertical.png',dpi=dpi); close



#RGR vs D
#D: sizez (in micrometers)
#spmcz=aggmassz+lpmz
#G: gd
#RGR=(alpha*fc*G*lpm*0.001*(lpm+aggmass/(1-porosity))/pho-beta*G^3/2*D^2)/aggmass
RGRz=(40*0.5*gd*lpmz*0.001*(aggmassz/0.02)/2600-14000*gd**1.5*(sizez/1e6)**2)/aggmassz #lpmz+
fig,ax=plt.subplots(1,1)
ax.set_xlabel('ESD (micrometer)')
ax.set_ylabel('RGR')
ax.scatter(sizez,log(RGRz))

savefig('logRGR.png'); close

#plt.show()







