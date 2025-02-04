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

ncfile='compare/coagulationbreakup.nc'#without differential settling
ncfile1='compare/diffset.nc'#'Results/v4h/results_vcurrent.nc' #sys.argv[1] #default output file
ncfile2='compare/sinking.nc' #'Results/vwithslackwater/results_vcurrent.nc'#'balemethod2.nc'

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
lpm1 = squeeze(ncv1['spm_spm'][:,lay])
lpm2 = squeeze(ncv2['spm_spm'][:,lay])
#agglpm = squeeze(ncv['agg_agglpm'][:,lay])
#totlpm=lpm+agglpm


#size layer SPMC
lpml = squeeze(ncv['spm_spm'][:,laysize]) 
lpml1 = squeeze(ncv1['spm_spm'][:,laysize]) #attention to ncv1
lpml2 = squeeze(ncv2['spm_spm'][:,laysize]) #attention to ncv2
#coagulationbreakup
lpmlmid = squeeze(ncv['spm_spm'][:,12]) #other layer in between
lpmlup = squeeze(ncv['spm_spm'][:,18]) #other layer in upper layer
lpmlbt = squeeze(ncv['spm_spm'][:,3]) #other layer in bottom layer
#+diffsett
lpmlmid1 = squeeze(ncv1['spm_spm'][:,12]) #other layer in between
lpmlup1 = squeeze(ncv1['spm_spm'][:,18]) #other layer in upper layer
lpmlbt1 = squeeze(ncv1['spm_spm'][:,3]) #other layer in bottom layer
#+sinking
lpmlmid2 = squeeze(ncv2['spm_spm'][:,12]) #other layer in between
lpmlup2 = squeeze(ncv2['spm_spm'][:,18]) #other layer in upper layer
lpmlbt2 = squeeze(ncv2['spm_spm'][:,3]) #other layer in bottom layer


h=3600./60 #/120. output was per 120s
tstart=int(2*h)#int(6.*h) #inital time for plotting
tend=tstart+int(39*h)#tstart+int(8*h)
tslice=slice(tstart,tend)
ptime=secs[tslice]/3600.


dpi=200
figure(figsize=(9,9*0.75))
lw=2.0
plot(ptime,lpm[tslice],'k',lw=lw+1,label='simulated SPM' )#1xTEP'total lit3''TSM')
plot(ptime,lpm1[tslice],'k--',lw=lw+1,label='+diffsett' )#+ differential settling
plot(ptime,lpm2[tslice],'k:',lw=lw+1,label='+sinking' )#+ sinking

plot(ptime,lpml[tslice],'r',lw=lw+1,label= 'size layer SPMC')
plot(ptime,lpml1[tslice],'r--',lw=lw+1,label= 'size layer + diffsett')
plot(ptime,lpml2[tslice],'r:',lw=lw+1,label= 'size layer + sinking')


#plot(ptime,lpmlmid[tslice],'pink',lw=lw+1,label= 'mid layer SPMC')#'with differential settling, mid layer
#plot(ptime,lpmlup[tslice],'orange',lw=lw+1,label= 'upper layer SPMC')#'with differential settling, mid layer
plot(ptime,lpmlbt[tslice],'purple',lw=lw+1,label= 'lower layer SPMC')#'mid layer

#plot(ptime,lpmlup1[tslice],'orange',linestyle='--',lw=lw+1,label= 'upper layer SPMC')#'with differential settling, mid layer
plot(ptime,lpmlbt1[tslice],'purple',linestyle='--',label= 'lower layer + diffset')#'with differential settling, mid layer
plot(ptime,lpmlbt2[tslice],'purple',linestyle=':',label= 'lower layer + sinking')#'with differential settling, mid layer

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
lpm=np.array(lpm)
lpml=np.array(lpml)
#print (min(totlpm[200:]),min(totlpml[200:]))



##vertical profile, SPMC g/m-3, log scale
dpi=400
##construct array for contour plot
lpmz = squeeze(ncv['spm_spm'][:,:])
z=squeeze(ncv['z'][:,:])
lpmz = transpose(lpmz)
timez=np.c_[secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs,secs]
timez=transpose(timez)
datatime=np.c_[data[:,0]-5,data[:,0]-5,data[:,0]-5]
dataz=np.c_[-(0.28-data[:,1]*0.01),-(0.28-data[:,1]*0.01),-(0.28-data[:,1]*0.01)] #from cm to m #-(0.28-23*0.01)
dataspmc=np.c_[data[:,2],data[:,2],data[:,2]]
fig,ax=plt.subplots(1,1)
cp = ax.contourf(timez/3600,transpose(z),np.log10(lpmz),cmap='Reds') #np.log10(spmcz) #norm=matplotlib.colors.LogNorm()(vmin=spmcz.min(), vmax=spmcz.max()))
#norm=matplotlib.colors.PowerNorm(gamma=1. / 2.),
#                       cmap='PuBu_r')
colorbar=fig.colorbar(cp, ax=ax, extend='max')
colorbar.set_label('SPMC',labelpad=-40, y=1.05, rotation=0) #log(SPMC)


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
gd = squeeze(ncv['agg_G'][:,:])
gd = transpose(gd)
fig,ax=plt.subplots(1,1)
ax.set_xlabel('time (h)')
ax.set_ylabel('water depth (m)')
cpg = ax.contourf(timez/3600,zesd,log(gd),cmap=cm.Blues)##,locator=ticker.LogLocator())
#levels=[0,0.001,0.01,0.1,1,5,10]
#norm = simple_norm(image, 'sqrt')
clbg=fig.colorbar(cpg, ax=ax)#,origin='lower',norm=norm)#,ticks=levels)
clbg.set_label('ln(G)',labelpad=-40, y=1.05, rotation=0)
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
savefig('gprofile.png',dpi=dpi); close
#print (len(gd),len(zesd1))
#print (np.max(gd))


#Dsize profile at certain time
fig,ax12=plt.subplots(1,1)
dpro = squeeze(ncv['agg_Dsize'][:,:]) #Dsize profile
#dpro = transpose(dpro)
ax12.plot(dpro[60*2]*1e6,zesd1[0],label='2 h')
ax12.plot(dpro[60*3]*1e6,zesd1[0],label='3 h')
ax12.plot(dpro[60*4]*1e6,zesd1[0],label='4 h')
ax12.set_xlabel('D ($\mu$m)')
ax12.set_ylabel('water depth (m)')
ax12.legend(loc='upper right')
savefig('Dsizeprofile.png',dpi=dpi); close
#print (len(gd),len(zesd1))
#print (np.max(gd))


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


G = squeeze(ncv['agg_G'][:,lay])
ws = squeeze(ncv['agg_ws'][:,lay])
ws1 = squeeze(ncv1['agg_ws'][:,lay])
ws2= squeeze(ncv2['agg_ws'][:,lay])

#changed to size layer
lpml = squeeze(ncv['spm_spm'][:,laysize])
Gl = squeeze(ncv['agg_G'][:,laysize])
wsl = squeeze(ncv['agg_ws'][:,laysize])
wsl1 = squeeze(ncv1['agg_ws'][:,laysize])
wsl2 = squeeze(ncv2['agg_ws'][:,laysize])

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
lpm1 = squeeze(ncv1['spm_spm'][:,laysize])
lpm2 = squeeze(ncv2['spm_spm'][:,laysize])
G = squeeze(ncv['agg_G'][:,laysize])
Gl = squeeze(ncv['agg_G'][:,laysize])
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
dsize=squeeze(ncv['agg_Dsize'][:,laysize])
dsize1=squeeze(ncv1['agg_Dsize'][:,laysize])
dsize2=squeeze(ncv2['agg_Dsize'][:,laysize])
dsizel=squeeze(ncv['agg_Dsize'][:,lay]) #Dsize for SPM layer
dsizel1=squeeze(ncv1['agg_Dsize'][:,lay]) #+differential settling
dsizel2=squeeze(ncv2['agg_Dsize'][:,lay]) #+sinking

##subplots 2 ssc and size
#fig, (ax1, ax3) = plt.subplots(2, sharex=True)
fig,ax1=plt.subplots(1,1)
ax1.set_ylabel('size $\mu$m')#, color=color

ax1.plot(ptime,1e6*dsize[tslice],label='1D dynamical $D$',color='red',lw=lw+1)
ax1.plot(ptime,1e6*dsize1[tslice],'r--',label='with differential settling',lw=lw+1)
ax1.plot(ptime,1e6*dsize2[tslice],'r:',label='with sinking',lw=lw+1)
#ax1.plot(ptime,1e6*dsizel[tslice],label='1D dynamical SPM layer $D$',color='blue',lw=lw+1)
#ax1.plot(ptime,1e6*dsizel1[tslice],'b--',label='with differential settling',lw=lw+1)
#ax1.plot(ptime,1e6*dsizel2[tslice],'b:',label='with sinking',lw=lw+1)
ax1.set_ylim()
ax1.legend(loc='upper right')

#ax3.plot(ptime,lpm[tslice],'k--',lw=lw+1,label='1D diagnostic $D$')#1xTEP'total lit3''TSM')
#ax3.plot(data[:,0]-5,data[:,2],'k^',ms=5.,label='data')#'TSM\nBale.ea2002')
#ax3.set_xlabel('time (h)')
#ax3.set_ylabel('SPM concentration (g m$^{-1}$)')
#ylim(0,4500)
#xlim(6,14)
#xlim(0,40)

s1datatime=sizedata[:,0]-5
s1datasize=sizedata[:,1] #data is in micrometer
#ax1.scatter(s1datatime, s1datasize,color='green',label='datasize')

ax1.set_ylim(20,)
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
coagulationlpm1=squeeze(ncv1['agg_coagulationlpm'][:,laysize]) #+ diffset

coagulationcalculated=14*0.5*Gl[tslice]*lpml[tslice]*1e-3*(dsizel[tslice])**2
###coagulation*1.d-3*lpm*Dsize**(4-self%fractal_dimension)
breakuplpm=squeeze(ncv['agg_Breakup'][:,laysize])
breakuplpm1=squeeze(ncv1['agg_Breakup'][:,laysize]) #+ diffset

breakupcalculated=13000*(Gl[tslice])**1.5*(dsizel[tslice]-1e-5)*(dsizel[tslice])**2

diffsetlpm=squeeze(ncv1['agg_diffsetlpm'][:,laysize])
diffsetlpmbt=squeeze(ncv1['agg_diffsetlpm'][:,1])
diffsetlpmtop=squeeze(ncv1['agg_diffsetlpm'][:,18])
diffsetlpmbt1=squeeze(ncv1['agg_diffsetlpm'][:,2])
###-self%breakup_factor*G**1.5d0*(Dsize-self%min_size)*Dsize**2 
#ax.plot(ptime,1e-3*agglpm[tslice]*breakuplpm[tslice],'--',lw=lw+1,label='breakup rate') #breakup rate in kg/m-3/s
#ax.plot(ptime,1e-3*coagulationlpm[tslice],'-,',lw=lw+1,label='coagulation rate') #ptime G  #coagulation rate in kg/m-3/s
#ax.plot(ptime,-breakupcalculated,'--',lw=lw+1,label='breakup rate') #breakuplpm[tslice] , #breakup rate dynamical size change  #ptime  #Gl is laysize G
#ax.plot(ptime,coagulationcalculated,'-,',lw=lw+1,label='coagulation rate') #coagulationlpm[tslice] coagulation rate dynamical size change

#ax.plot(ptime,-breakuplpm[tslice],'--',lw=lw+1,label='breakup rate diag') #Gl[tslice]**1.5
#ax.plot(ptime,coagulationlpm[tslice],'-,',lw=lw+1,label='coagulation rate diag') #Gl[tslice] #totlpm1[tslice]
ax.plot(ptime,diffsetlpm[tslice],'-,',lw=lw+1,label='diffset rate diag') #Gl[tslice] #totlpm1[tslice]
ax.plot(ptime,diffsetlpmbt[tslice],'-,',lw=lw+1,label='diffset rate bottom') #Gl[tslice] #totlpm1[tslice]
ax.plot(ptime,diffsetlpmtop[tslice],'-,',lw=lw+1,label='diffset rate top') #Gl[tslice] #totlpm1[tslice]
ax.plot(ptime,diffsetlpmbt1[tslice],'-,',lw=lw+1,label='diffset rate near bottom') #Gl[tslice] #totlpm1[tslice]

legend()
savefig('coagulationbreakup.png',dpi=dpi); close
#print (breakuplpm[tslice]/coagulationlpm[tslice]) #*(Gl[tslice])**0.5
#print (breakupcalculated/coagulationcalculated)

f=figure(figsize=(9,4.5),dpi=dpi)
#f.subplots_adjust(bottom=0.15)
lw=2.0
#plot(ptime,size(lpm[tslice],agglpm[tslice]),'k-',lw=lw)
plot(ptime,1000*ws[tslice],'k-',lw=lw, label='size layer')
plot(ptime,1000*ws1[tslice],'k--',lw=lw, label='size layer + diffsett')
plot(ptime,1000*ws2[tslice],'k:',lw=lw, label='size layer + sinking')
#plot(ptime,1000*wsl[tslice],'b',lw=lw, label='SPMC layer')
#plot(ptime,1000*wsl1[tslice],'b--',lw=lw, label='SPMC layer + diffset')
xlabel('time [h]')
ylabel(r'sinking velocity [mm/s]')
#ylim(-20,0)
legend()
savefig('Ws.png',dpi=dpi); close()
ws=np.array(ws)
wsl=np.array(wsl)
#print (max(ws[500:]),max(wsl[500:]))


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

#simple comparison
test = squeeze(ncv['agg_Dsize'][240:,laysize])#laysize
test1 = squeeze(ncv1['agg_Dsize'][240:,laysize])
test2 = squeeze(ncv2['agg_Dsize'][240:,laysize])
test3=test1-test
test4=test2-test
#print (test.shape)
print ("max D:",np.max(test)*1e6,"min D:",np.min(test)*1e6)
print ("max D + diffsett:",np.max(test1)*1e6,"min D + diffsett:",np.min(test1)*1e6, "max%:",(np.max(test1)-np.max(test))*100/np.max(test), "min%:",(np.min(test1)-np.min(test))*100/np.min(test))
print ("max D + sinking:",np.max(test2)*1e6,"min D + sinking:",np.min(test2)*1e6,"max%:",(np.max(test2)-np.max(test))*100/np.max(test),"min%:",(np.max(test2)-np.max(test))*100/np.min(test))

print ("+/- diffsett max:",np.max(test3)*1e6)
print ("+/- sinking max:",np.max(abs(test4))*1e6)











