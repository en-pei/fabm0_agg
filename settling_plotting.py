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
from warnings import filterwarnings
from scipy.signal import argrelextrema
import seaborn as sns
import pandas as pd
from scipy import stats
import statistics
from scipy.optimize import curve_fit
import matplotlib.patches as patches
import matplotlib.lines as mlines
#import itertools
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import FunctionTransformer
from sklearn import linear_model
import warnings
from datetime import timedelta
import string
import matplotlib.patheffects as PathEffects
import matplotlib.colors as mcolors
from matplotlib.patches import Patch



warnings.simplefilter('ignore',DeprecationWarning)



##define dataset and model result directory
datadir = '/gpfs/home/lie/setups/bale2002/Plot/data' #directory for data 
folder='/gpfs/work/lie/run/par/plotting/breakup_factor-1310-coagulation_rate-1.583-kws-347/' # directory for model results 
baleref='/gpfs/work/lie/run/par/plotting/breakup_factor-1310-coagulation_rate-1.583-kws-347/bale.nc'  #lab base run folder
fettref='/gpfs/work/lie/run/par/6oszmay/testplot/bale2002_.nc' #field base run folder 


##### get model run nc file results 
def process_data(filename: str): #, lay: int, laysize: int, layflow: int
    lay=16
    laysize=10
    layflow=9
    ncv = netCDF4.Dataset(filename).variables
    utime = cftime.utime(ncv['time'].units)
    time = utime.num2date(ncv['time'][:])
    secs = ncv['time'][:]
    lpm = np.squeeze(ncv['spm_spm'][:, lay])
    dsize = np.squeeze(ncv['agg_Dsize'][:, laysize])
    esd = dsize  # You may adjust this if needed
    G = np.squeeze(ncv['agg_G'][:, lay])
    uu = np.abs(np.squeeze(ncv['u'][:, layflow])) * np.sqrt(2)
    return secs, lpm, esd,G, uu



secs0, lpm0, esd0, G0, uu0 = process_data(folder+'bale.nc') #base run
secs1, lpm1, esd1, G1, uu1 = process_data(folder+'temp/0kb.nc') #
secs2, lpm2, esd2, G2, uu2 = process_data(folder+'baleepsexp.nc') #with degradation
secs3, lpm3, esd3, G3, uu3 = process_data(folder+'temp/0ks.nc') #




n=2 #plot per n point in data 
mks=70 #20 #markersize for mobaxterm
plt.rcParams.update({'font.size': 20})  # 22
plt.rcParams.update({'lines.linewidth':3}) #3
###all of the time of this data needed /20*16+4 #got the wrong interval first 
bale3500spmc=asarray(numpy.loadtxt('%s/spmc_bale3500_new.txt'%datadir,delimiter=','))  #','
bale3500spmcf=asarray(numpy.loadtxt('%s/flow_balespmc3500_new.txt'%datadir,delimiter=','))
bale3500d50=asarray(numpy.loadtxt('%s/d50_bale3500_new.txt'%datadir,delimiter=','))
bale3500d50f=asarray(numpy.loadtxt('%s/flow_baled503500_new.txt'%datadir,delimiter=','))#newly datathieved 20221004

oszspmmay=asarray(numpy.loadtxt('%s/oszspmmay.txt'%datadir,delimiter=','))
oszd50may=asarray(numpy.loadtxt('%s/oszd50may.txt'%datadir,delimiter=','))
oszflowmay=asarray(numpy.loadtxt('%s/oszflowmay.txt'%datadir,delimiter=','))


tsmdatat=bale3500spmc[:,0]/20*16+4 -0.95-0.1-0.6 +0.2
d50datat=bale3500d50[:,0]/20*16+4 -0.95-0.1-0.6 +0.2 
flowdatat=bale3500spmcf[:,0]/20*16+4 -0.95-0.1-0.6 +0.2


tsmdata=bale3500spmc[0:,1] 
d50data=bale3500d50[0:,1] 
flowdata=bale3500spmcf[0:,1] 


tsmdata1=oszspmmay[:,1] 
d50data1=oszd50may[:,1]
flowdata1=oszflowmay[:,1] 


##Get maerz&wirtz 2009 data from fig 4.
maerz09datadfeb=asarray(numpy.loadtxt('%s/maerz2009_d50_feb.txt'%datadir,delimiter=','))
maerz09dt=maerz09datadfeb[:,0]*8+8
maerz09d=maerz09datadfeb[:,1]
colorm='r'

data = asarray(numpy.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=','))
sizedata = asarray(numpy.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=','))
s1datatime=sizedata[:,0]-5#-0.25
s1datasize=sizedata[:,1] #data is in micrometer
datatime=data[:,0]-5.#-0.25
dataspmc=data[:,2]



######################### fig4 0kb and 0ks ################################
fig, axs = plt.subplots(2,1,sharex=True)
axs[0].plot(secs0/3600,esd0*1e6,color='k',lw=3,label='settling+breakup')#,label='base run') #'
axs[0].plot(secs1/3600,esd1*1e6,color='orange',lw=2.5,label='only settling')#,label='only settling')  #'b'
axs[0].plot(secs3/3600,esd3*1e6,color='r',alpha=0.8,lw=2,label='only breakup')#,label='only breakup')  #'b'
axs[1].plot(secs0/3600,lpm0,color='k',lw=3) #base run
axs[1].plot(secs1/3600,lpm1,color='orange',lw=2.5)
axs[1].plot(secs3/3600,lpm3,color='r',alpha=0.8,lw=2)
axs[0].set_ylabel('mean size (μm)') #
axs[0].set_xlim(7.1,16)
axs[0].set_ylim(50,230)
axs[1].set_ylim(0,)
 #7.1,16 7,12   7,15    5,13
axs[1].set_ylabel('SPMC (mg/L)' ) #g m$^{-3}$  

###twin axis for showing flow velocity
coloru='grey' #'orange' #color of uu. flow velocity
ax5 = axs[0].twinx()
ax5.plot(secs0/3600,uu0,linewidth=0.8,color=coloru,alpha=0.5,label='tide') #try twin yaxis '--',
ax5.set_ylim(0,0.65) #0,0.75
ax5.set_ylabel('velocity (m/s)',fontsize=20 )#('ESD ($\mu$m)') velocity  velocity (m s$^{-1}$)
ax5.title.set_color(coloru)
ax5.yaxis.label.set_color(coloru)
ax5.tick_params(axis='y', colors=coloru)

axs[1].set_xlabel('time (h)')

#axs[1].scatter(tsmdatat[::n],tsmdata[::n], s=mks, facecolors='none', edgecolors='k')
#axs[0].scatter(d50datat[::n]-4,d50data[::n], s=mks, facecolors='none', edgecolors='k',label='observations') #Bale data 'measurement'
axs[0].legend(loc='upper center',bbox_to_anchor=(0.5,1.4),framealpha=0.9,ncol=4)
#axs[1].legend(loc='upper center',bbox_to_anchor=(1.3,0.8),framealpha=0.9,ncol=1)

####### end of plotting only settling and only breakup ######
plt.show()
exit()


#################### box plot  #############################
df=pd.DataFrame({'time': d50datat[5:]-4,'size': d50data[5:]})
df['time']=pd.to_datetime(df['time']*3600,unit='s',origin='unix')
#df['mark']='data'
df['mark']='observations'
#df['4h']=pd.to_datetime(df['time']).dt.floor('4H')

dfmodel=pd.DataFrame({'time': secs2[0:60*16],'size': esd2[0:60*16]*1e6}) #decomposition model run
dfmodel['time']=pd.to_datetime(dfmodel['time'],unit='s',origin='unix')
#dfmodel['mark']='model'
dfmodel['mark']='model with decomposition'

##df1=pd.merge(df,dfmodel,how='left',on='time')
df1=pd.concat([df,dfmodel],ignore_index=True)
df1=df1.sort_values(by='time',ascending=True)
#sns.catplot(data=df.assign(group=df['time'].dt.floor('4h').dt.time),x='group',y='size',kind='box')

fig, axs = plt.subplots() #2,1,sharex=True,sharey=True
dfnew=df1.assign(group=df1['time'].dt.floor('4h').dt.time)
dfnew['group1']=(pd.to_datetime(dfnew['group'],format='%H:%M:%S')+pd.DateOffset(hours=2)).dt.hour
#dfnew['group1']=dfnew['group'].apply(pd.Timestamp)+pd.DateOffset(hours=2)
#print (dfnew['group1'])

bx=sns.boxplot(data=dfnew,x='group1',y='size',hue='mark',ax=axs) #group
#bx=sns.boxplot(data=df1.assign(group=df1['time'].dt.floor('4h').dt.time),x='group',y='size',hue='mark',ax=axs)
sns.move_legend(bx,"upper left",bbox_to_anchor=(1,1))
#sns.boxplot(data=df.assign(group=df['time'].dt.floor('4h').dt.time),x='group',y='size',palette='Greens',ax=axs[0])
#sns.boxplot(data=dfmodel.assign(group=dfmodel['time'].dt.floor('4h').dt.time),x='group',y='size',palette='Blues',ax=axs[1])
#plt.show()
#exit()

df.set_index('time',inplace=True)
meandata=df.resample('4H',closed='right').mean()
meandata.reset_index(inplace=True)
meandata['hour']=meandata['time'].dt.hour

maxdata=df.resample('4H',closed='right').max()
maxdata.reset_index(inplace=True)
maxdata['hour']=maxdata['time'].dt.hour



#######################Fig 6 size drifting ########
### to make a box pt for data and model comparisions###
#################### fig4 and fig6 plot code ####### decomposition #########
plt.close('all')
fig, axs = plt.subplots(2,1,sharex=False)  ####something funky with the timeline here when share x axis
#axs[0].axhline(y=130,color='k',ls='--',lw=0.8) #reference line
axs[0].plot(secs0/3600,esd0*1e6,color='k',lw=2,label='base run') #'secs0/3600 pd.to_datetime(secs0*6,unit='s',origin='unix')
axs[0].plot(secs2/3600,esd2*1e6,color='mediumseagreen',lw=3,label='with decomposition')  #'b'
#axs[0].fill_between(secs2/3600,esd2*1e6,esd0*1e6,color='none',edgecolor='k',hatch='\\\\') #secs2/3600 color='grey',alpha=0.2,#shaded area
axs[0].set_ylabel('simulated size (μm)') #μ
##axs[1].legend(loc='upper right',bbox_to_anchor=(1.25,0.8),framealpha=0.9)
axs[0].xaxis.set_major_locator(ticker.MultipleLocator(4))
axs[0].set_xticklabels([])
axs[0].legend(loc="upper center",bbox_to_anchor=(0.5,1.2),ncol=2)
axs[0].set_xlim(0,16)
axs[0].set_ylim(87,215)
#axs[0].set_xlim(0.1,23.9)
axs[1].axhline(y=130,color='gray',ls='--',lw=0.8,zorder=0) #reference line
dfnew['hour']=dfnew['group1']
palette=['mediumaquamarine','lightgrey'] #'mediumseagreen','lightgrey'

#bx=sns.boxplot(data=dfnew,x='hour',y='size',hue='mark',ax=axs[1],palette=palette) #group #original works fine 20241222

bx = sns.boxplot(
    data=dfnew,
    x='hour',
    y='size',
    hue='mark',
    ax=axs[1],
    palette=palette,
    width=0.5,  # Make the boxes thinner (adjust the box width)
    medianprops={"linewidth": 5, "color": "black"},  # Thicker and customized median line
#    showmeans=True,
#    meanprops={"marker": "o", "markerfacecolor": "red", "markeredgecolor": "black", "markersize": 8}
)

# Apply alpha to the fill color of the boxes
alpha_value=0.3
for patch in bx.patches:
    facecolor = patch.get_facecolor()  # Get current facecolor
    patch.set_facecolor((*facecolor[:3], alpha_value))  # Add alpha to the RGB color

#bx=sns.boxplot(data=df1.assign(group=df1['time'].dt.floor('4h').dt.time),x='group',y='size',hue='mark',ax=axs[1])
axs[1].set_ylabel('size range during a cycle') #
axs[1].set_xlabel('time (h)') #
axs[1].set_ylim(87,172)

# Create updated legend handles with alpha-adjusted colors
handles, labels = bx.get_legend_handles_labels()  # Get current legend handles and labels
updated_handles = [
    Patch(facecolor=(*sns.color_palette(palette)[i][:3], alpha_value), edgecolor='black') 
    for i in range(len(palette))
]

# Update the legend with the new handles while keeping the exact location and style
sns.move_legend(
    bx,
    "upper center",
    bbox_to_anchor=(0.5, 1.2),
    ncol=2,
    title=None,
    handles=updated_handles,
    labels=labels
)

######################### end of decomposition plot fig.6 now##########################



################# read model lab #############################

lay=16 ##layer for SPM concentration
laysize=9 ##layer for ESD
layflow=9 #lay flow
ncv=netCDF4.Dataset(baleref).variables

utime=cftime.utime(ncv['time'].units)
time=utime.num2date(ncv['time'][:])
secs=ncv['time'][:]
lpm=squeeze(ncv['spm_spm'][:,lay])
#lpmsize = squeeze(ncv['spm_spm'][:,laysize])  #in size layer
dsize=squeeze(ncv['agg_Dsize'][:,laysize]) #  dynamical size
esd=dsize #squeeze(ncv['agg_esd'][:,laysize])  #dsize #dsize is for dynamical
G = squeeze(ncv['agg_G'][:,lay])
uu = abs(squeeze(ncv['u'][:,layflow]))*sqrt(2) #flow velocity


################## read model field ####################
#'/gpfs/work/lie/oszmay_newphysics_minErrRange/bale2002.nc'
#'/gpfs/work/lie/run/par/6oszmay/model/bale2002.nc'  
lay1=6  #layer for SPM concentration
laysize1=6  #layer for ESD
layflow1=6 #lay flow
ncv1=netCDF4.Dataset(fettref).variables

utime1=cftime.utime(ncv1['time'].units)
time1=utime1.num2date(ncv1['time'][:])
secs1=ncv1['time'][:]
lpm1=squeeze(ncv1['spm_spm'][:,lay1])
#lpmsize = squeeze(ncv1['spm_spm'][:,laysize1])  #in size layer
dsize1=squeeze(ncv1['agg_Dsize'][:,laysize1]) #  dynamical size
esd1=dsize1 #squeeze(ncv['agg_esd'][:,laysize1])  #dsize #dsize is for dynamical
G1= squeeze(ncv1['agg_G'][:,lay1])
uu1 = abs(squeeze(ncv1['u'][:,layflow1]))*sqrt(2) #flow velocity

################### plotting fig2#####################################
#################### 6 panel plots for validation and stackplot #####
plt.close('all')
fig, axs = plt.subplots(3,2,  sharex='col')  #4 (5, sharex=False)
#fig, axs = plt.subplot_mosaic("ab;cd;ef")  
#_=[axs[i].sharex(axs['a']) for i in 'ace']

for nn, ax in enumerate(axs.flat):
	ax.text(0.05,0.85,string.ascii_lowercase[nn],transform=ax.transAxes,size=20,weight='bold') #add subplot numbers 


axs[0,0].plot(secs/3600,esd*1e6,label='model',color='k') #'i.name',+0.27#in m #ESD t  #[0,0]
axs[0,0].scatter(d50datat[::n]-4,d50data[::n], s=mks, facecolors='none', edgecolors='k',label='observations') #Bale data 'measurement'

axs[0,1].plot(secs1/3600,esd1*1e6,color='k')
axs[0,1].scatter(d50datat1[::n],d50data1[::n], s=mks, facecolors='none', edgecolors='k',label='Bale data')
axs[0,1].set_ylim(0,300)

########### put SPMC in the lowest panel ###
axs[2,0].plot(secs/3600,lpm,color='k')#,label='model') # i.name ,color='firebrick', totlpm    #+0.27 #SPMC	vs t
axs[2,0].scatter(tsmdatat[::n],tsmdata[::n], s=mks, facecolors='none', edgecolors='k')

axs[2,1].plot(secs1/3600,lpm1,color='k')  #oszmay spmc
axs[2,1].scatter(tsmdatat1[::n],tsmdata1[::n], s=mks, facecolors='none', edgecolors='k')

colorm='dimgrey'
coloru='violet' #'orange' #color of uu. flow velocity

ax5 = axs[0,1].twinx()
ax5.plot(secs1/3600,uu1,linewidth=2,color=coloru,alpha=0.5,label='tide') #try twin yaxis,'--',
ax5.set_ylim(0,0.65) #0,0.75
ax5.set_ylabel('velocity (m/s)',fontsize=20 )#('ESD ($\mu$m)') velocity  velocity (m s$^{-1}$)
ax5.title.set_color(coloru)
ax5.yaxis.label.set_color(coloru)
ax5.tick_params(axis='y', colors=coloru)


ax4 = axs[0,0].twinx()
ax4.scatter(maerz09dt,maerz09d,s=mks, marker='^', color='gainsboro' ,edgecolor='dimgrey',label='observations (Wadden Sea)')  #'Wadden data' whitesmoke  Lunau (2006) Lunau data
############size data from bale2002#######################################################################
###########################ax4 maerz09 data##########################################
#ax4 = axs[1].twinx()
##fig.subplots_adjust(right=0.75)
#ax4.spines['right'].set_position(('axes', 1.2))
#ax4.set_frame_on(True)
#ax4.patch.set_visible(False)
#ax4.scatter(maerz09dt,maerz09d, marker='^', color=colorm,label='Lunau(2006)')
ax4.set_ylim(10,50)
#ax4.set_ylabel('Lunau $D_{50}$ ($\mu$m)',color=colorm) #'%s Thing' % color
ax4.tick_params(axis='y', colors=colorm)


ax2 = axs[0,0].twinx()
ax2.plot(secs/3600,uu,linewidth=2,color=coloru,alpha=0.5,label='tide') #try twin yaxis '--',
ax2.set_ylim(0,0.65) #0.55 0,0.75
#ax2.set_ylabel('velocity (m s$^{-1}$)',fontsize=20 )#('ESD ($\mu$m)') velocity 
#ax2.title.set_color('grey')
#ax2.yaxis.label.set_color('grey')
ax2.set_yticklabels([]) #remove the number
ax2.set_yticks([])
#ax2.tick_params(axis='y', colors='grey')


axs[0,0].set_ylim(0,300)
axs[0,0].set_xlim(7.1,16) #7,12   7,15    5,13
#axs[2,0].set_xlabel('time (h)')

axs[0,0].set_ylabel('size (μm)')
axs[0,0].xaxis.set_major_locator(plt.MaxNLocator(3)) #how many ticks for x axis

axs[1,0].set_ylabel('process \ncontribution' ) #g m$^{-3}$  , fontsize=20  'SPMC (% of max)' #size factors

axs[2,0].set_ylabel('SPMC (mg/L)' ) #g m$^{-3}$  , fontsize=20  'SPMC (% of max)'
axs[2,0].set_ylim(0,)
axs[2,1].set_ylim(0,50) #(0,120) #
axs[2,1].set_xlim(5,18)

#axs[1,1].set_xlim(5,18)

handles1, labels1 = axs[0,0].get_legend_handles_labels()
handles4, labels4 = ax4.get_legend_handles_labels()

# Combine handles and labels
handles = handles1 + handles4
labels = labels1 + labels4

# Create a legend without a box
#axs[1,0].legend(handles, labels, loc='upper right') #, frameon=False
#axs[0, 0].legend(handles, labels, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.3), bbox_transform=axs[0, 0].transAxes) #bbox_to_anchor=(1.1, 2.1)

#fig.tight_layout()



dDoszmay=pickle.load(open('/gpfs/work/lie/bale2002/Plot/result_to_download/dDi_oszmay.pkl','rb'))
dDbale=pickle.load(open('/gpfs/work/lie/bale2002/Plot/result_to_download/dDi_bale.pkl','rb'))


dDbale['time']=dDbale['time']+4
dDbale['sumd']=abs(dDbale.coagulation)+abs(dDbale.diffset)+abs(dDbale.breakup)+abs(dDbale.settling)+abs(dDbale.diffusion)

dDoszmay['sumd']=abs(dDoszmay.coagulation)+abs(dDoszmay.diffset)+abs(dDoszmay.breakup)+abs(dDoszmay.settling)+abs(dDoszmay.diffusion)
print (dDbale)
print (dDbale.info())
print (dDoszmay)

####### put stackplot in the middle panel
axs[1,0].stackplot(dDbale.time,dDbale.coagulation/dDbale.sumd,dDbale.diffset/dDbale.sumd,baseline ='zero',colors =['mediumseagreen', 'b']) 
axs[1,0].stackplot(dDbale.time,dDbale.settling/dDbale.sumd,-dDbale.breakup/dDbale.sumd,baseline ='zero',colors =['orange','r']) 


axs[1,1].stackplot(dDoszmay.time,dDoszmay.coagulation/dDoszmay.sumd,dDoszmay.diffset/dDoszmay.sumd,baseline ='zero',colors =['mediumseagreen', 'b']) 
axs[1,1].stackplot(dDoszmay.time,dDoszmay.settling/dDoszmay.sumd,-dDoszmay.breakup/dDoszmay.sumd,baseline ='zero',colors =['orange','r']) 

axs[1,0].set_ylim(-0.59,0.41) #-0.7,0.51 -0.75,0.75
axs[1,1].set_ylim(-0.59,0.41)

lwnr=8

axs[1,0].plot([], [],lw=lwnr, color ='mediumseagreen', #axs[1,0]
         label ='coagulation (turbulence)')  #'s',
axs[1,0].plot([], [], lw=lwnr,color ='b',
         label ='coagulation (settling)') #'differential settling'
axs[1,0].plot([], [], lw=lwnr,color ='orange',
         label ='preferential settling') #sinking
axs[1,0].plot([], [],lw=lwnr, color ='r',
         label ='breakup')
#axs[0,0].plot([], [], lw=lwnr,color ='lightgrey',
#         label ='diffusion')

axs[0,1].legend(handles, labels, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.8),bbox_transform=axs[0, 0].transAxes) #0.5,1.5 ##legend for model and data
axs[1,0].legend(loc='upper center',bbox_to_anchor=(0.5, 8.01),ncol=4) #### legend for stackplot 0.5, 6.01
#axs[1,0].legend(loc='center right',bbox_to_anchor=(6.5, 1),ncol=2) #### legend for stackplot



#for n,(key,ax) in enumerate(axs.items()):
#for ho in axes:
#	for i in ho:
#	ax.txt(-0.1,1.1,key,transform=ax.transAxies,size=20,weight='bold')
plt.subplots_adjust(wspace=0.35) #0.4  # Adjust the value as needed
plt.subplots_adjust(hspace=0.05) #0.15 
fig.text(0.5,0.04, 'time (h)',ha='center')
#plt.show()
#exit()
################### end of plotting 6 subplots ###


refdir='/gpfs/work/lie/p1ref'
############ Fig 3. Sensitivity ########################
############ make sensitivity plot of different Amp ####################
def flatten(l):
    return [item for sublist in l for item in sublist]

plt.close('all')
#fig, axs = plt.subplots(1,2)  #2,2,  sharex='col' 4 (5, sharex=False)




fr=75*6 #59*6 #90*6 #
to=141*6 #122*6 #90*6  #
lfr=8 #5 #8 #5 #0
lto=10 #7 #10 #7 #-1
#ncvi=ncv #only one file 

filelist=[]
gs=[]
ratios=[]
markers=[]
us=[]
layf=9 #6 #14 #6





for i in os.scandir(refdir+'/glab'): #/glab gfield
	i.is_dir()
	ncfilei=i.path#+'*.nc'  #+'*.nc' #/bale2002_.nc '/smallernf.nc' # #'/largerkc.nc' #'/bale2002.nc'  ##!!!!!!!!!!!!!!!###
	base_name = os.path.splitext(i.name)[0]
	filelist.append(base_name) #i.name
	#print (base_name)

	ncvi=netCDF4.Dataset(ncfilei).variables
	ratiocali=-squeeze(ncvi['agg_sinkinglpm'][fr:to,lfr:lto])/squeeze(ncvi['agg_Breakup'][fr:to,lfr:lto]) #fr:to

	ratioi=np.mean(ratiocali,axis=1)#[:-1]
	ratioi=list(ratioi)
	ratios.append(ratioi) 

	G4h = squeeze(ncvi['agg_G'][fr:to,lfr:lto])#.astype(bool)  #
	gmean=np.mean(G4h, axis=1)  #g4hall   max
	gcycle=list(gmean)
	gs.append(gcycle) #   gcycle*cyclen          #*5   5*4h cycle

	markeri=[base_name]*len(G4h) #ratiocali  ##Amplitude of M2 
	markers.append(markeri)
	uu = abs(squeeze(ncvi['u'][fr:to,layf]))*sqrt(2) #flow velocity 6 for field
	us.append(uu)

print (filelist)

ratios=flatten(ratios)
markers=flatten(markers)
gs=flatten(gs)
us=flatten(us)

#dmeans=flatten(dmeans)
#cmeans=flatten(cmeans)

#ratiocbs=flatten(ratiocbs)
#depths=flatten(depths)

dict = {'G (s$^{-1}$)': gs, 'ratio':ratios,'Amplitude(M2)':markers,'velocity':us}    #'size': dmeans,'SPMC': cmeans, 

print (len(gs),len(ratios),len(markers))
df = pd.DataFrame(dict)
df=df.sort_values(by='Amplitude(M2)',ascending=False)
#print (df)


def sensitivity3(filename,layf,fr,to):#9,6 
	#fr=75*6 #59*6 #90*6 #
	#to=141*6 #122*6 #90*6  #

#ncvi=ncv #only one file 
	filelist=[]
	gs=[]
	ds=[]
	ratios=[]
	markers=[]
	us=[]
	ts=[]
	layflow=layf #6 #14 #6
	lfr=layflow-1#8 #5 #8 #5 #0
	lto=layflow+1 #10 #7 #10 #7 #-1
	refdir='/gpfs/work/lie/p1ref'
	for i in os.scandir(refdir+filename): #'/glab' #/glab gfield
		i.is_dir()
		ncfilei=i.path#+'*.nc'  #+'*.nc' #/bale2002_.nc '/smallernf.nc' # #'/largerkc.nc' #'/bale2002.nc'  ##!!!!!!!!!!!!!!!###
		base_name = os.path.splitext(i.name)[0]
		filelist.append(base_name) #i.name
		#print (base_name)

		ncvi=netCDF4.Dataset(ncfilei).variables
		ratiocali=-squeeze(ncvi['agg_sinkinglpm'][fr:to,lfr:lto])/squeeze(ncvi['agg_Breakup'][fr:to,lfr:lto]) #fr:to

		ratioi=np.mean(ratiocali,axis=1)#[:-1]
		ratioi=list(ratioi)
		ratios.append(ratioi) 
	
		G4h = squeeze(ncvi['agg_G'][fr:to,lfr:lto])#.astype(bool)  #
		gmean=np.mean(G4h, axis=1)  #g4hall   max
		gcycle=list(gmean)
		gs.append(gcycle) #   gcycle*cyclen          #*5   5*4h cycle

		D4h = squeeze(ncvi['agg_Dsize'][fr:to,lfr:lto]*1e6)#.astype(bool)  #
		dmean=np.mean(D4h, axis=1)  
		dcycle=list(dmean)
		ds.append(dcycle) 
	
		markeri=[base_name]*len(G4h) #ratiocali  ##Amplitude of M2 
		markers.append(markeri)
		uu = abs(squeeze(ncvi['u'][fr:to,layflow]))*sqrt(2) #flow velocity 6 for field
		us.append(uu)


		t=ncvi['time'][fr:to]
		ts.append(t)
	ratios=flatten(ratios)
	markers=flatten(markers)
	gs=flatten(gs)
	ds=flatten(ds)
	us=flatten(us)
	ts=flatten(ts)
	dict = {'G (s$^{-1}$)': gs, 'ratio':ratios,'Amplitude(M2)':markers,'velocity':us, 'size':ds,'time':ts}    #'size': dmeans,'SPMC': cmeans, 
	#print (len(gs),len(ratios),len(markers))
	df = pd.DataFrame(dict)
#	df=df.sort_values(by='time',ascending=False) #False
	df=df.sort_values(by='Amplitude(M2)',ascending=False) #False
	return df
print ('testing function')
#dflab=sensitivity3('/glab',9,75*6,141*6)
dflab=sensitivity3('/glab',9,75*6,141*6)
dffield=sensitivity3('/gfield',6,59*6,122*6)

dflabkc=sensitivity3('/kc',9,75*6,120*6)
dflabks=sensitivity3('/ks_new',9,75*6,121*6) 
custom_palette={'gold','c','k'}

plt.close('all')
#fig, ax = plt.subplots() #plot for kc sensitivty
fig, ax = plt.subplots(1,2,sharey=True,sharex=True)  #2,2,  sharex='col' 4 (5, sharex=False)
#sns.lineplot(data=dflabkc, x='G (s$^{-1}$)',y='size',hue='Amplitude(M2)',alpha=1,palette=custom_palette,ax=ax) #,legend=False marker='o',edgecolor='none',
#sns.scatterplot(data=dflabkc, x='G (s$^{-1}$)',y='size',hue='Amplitude(M2)',alpha=1,marker='o',edgecolor='none',s=80,palette=custom_palette,ax=ax) #,legend=False 
sns.scatterplot(data=dflabks, x='G (s$^{-1}$)',y='size',hue='Amplitude(M2)',alpha=1,marker='o',edgecolor='none',s=50,palette=custom_palette,ax=ax[0]) #,legend=False 
sns.scatterplot(data=dflabkc, x='G (s$^{-1}$)',y='size',hue='Amplitude(M2)',alpha=1,marker='o',edgecolor='none',s=50,palette=custom_palette,ax=ax[1]) #,legend=False 
ax[0].legend(title='settling parameter $k_s$',labels=['+','base run','-'],ncol=3,loc='upper center',bbox_to_anchor=(0.5,1.3),fancybox=True) #'2x','1x','0.5x'labels=['2x','1x','0.5x'],
ax[1].legend(title='coagulation parameter $k_c$',ncol=3,labels=['+','base run','-'],loc='upper center',bbox_to_anchor=(0.5,1.3),fancybox=True) #'2x','1x','0.5x'labels=['2x','1x','0.5x'],
ax[0].set_xlabel('turbulence $G$ (s$^{-1}$)')
ax[0].set_ylabel('mean size ($\mu$m)')
#plt.show()
#exit()


xplot='velocity' #'G (s$^{-1}$)' #'G (s$^{-1}$)' #'time (h)' #
yplot='ratio' #'size' #'SPMC' #'SPMC'

plt.close('all')
fig, ax = plt.subplots(1,2,sharey=True,sharex=True)  #2,2,  sharex='col' 4 (5, sharex=False)
sns.lineplot(
    data=dflab, x=xplot,y=yplot,hue='Amplitude(M2)',lw=5,alpha=0.7,palette=custom_palette,ax=ax[0],legend=False)

sns.lineplot(
    data=dffield, x=xplot,y=yplot,hue='Amplitude(M2)',lw=5,alpha=0.7,palette=custom_palette,ax=ax[1],legend=False)

ax[0].text(0.1,4e3,'lab,h=0.28m',size=20) 
ax[1].text(0.1,4e3,'field,h=20m',size=20) 

ax[1].legend(title='relative tidal amplitude',labels=['2x','1x','0.5x'],frameon=False)
ax[0].set_yscale('log')
ax[0].set_xlabel('flow velocity (m/s)')
ax[1].set_xlabel('velocity (m/s)')

ax[0].set_ylabel('settling-to-breakup ratio')
ax[0].axhline(y=1, color=".7", dashes=(2, 1),lw=0.5, zorder=0)
ax[1].axhline(y=1, color=".7", dashes=(2, 1),lw=0.5, zorder=0)
ax[1].set_yscale('log')
#plt.show()
#exit()
############# end of plotting 2 case 3 sensitivity plot and function ####################




## ampitude for lab #####
df['Amplitude(M2)'].replace(['0.075'],'0.5x',regex=True,inplace=True) # '0.075'
df['Amplitude(M2)'].replace(['0.15'],'1x',regex=True,inplace=True) #'0.15'
df['Amplitude(M2)'].replace(['0.30'],'2x',regex=True,inplace=True) #'0.30'

#### amplitude for field ###
#df['Amplitude(M2)'].replace(['0.25'],'0.5x',regex=True,inplace=True) # '0.075'
#df['Amplitude(M2)'].replace(['0.5'],'1x',regex=True,inplace=True) #'0.15'
#df['Amplitude(M2)'].replace(['1.0'],'2x',regex=True,inplace=True) #'0.30'
#df['Amplitude(M2)'].replace(['1xx'],'0.5x',regex=True,inplace=True) # somehow some buggy turn 0.5

df=df.rename(columns={'Amplitude(M2)':'Amplitude ref'})


print (df['Amplitude ref'].unique())
print (df.info())

custom_palette={'0.5x':'pink','1x':'red','2x':'maroon'}
custom_palette1={'0.5x':'lightsteelblue','1x':'blue','2x':'navy'}
fig, ax = plt.subplots(2,1,sharex=True)  #2,2,  sharex='col' 4 (5, sharex=False)
sns.lineplot(
    data=df, x=xplot,y=yplot,hue='Amplitude ref',lw=5,alpha=0.7,palette=custom_palette,ax=ax[0],legend=False)
ax[0].set_yscale('log')
ax[0].axhline(y=1, color=".7", dashes=(2, 1),lw=0.5, zorder=0)
sns.lineplot(
    data=df, x=xplot,y='G (s$^{-1}$)',hue='Amplitude ref',lw=5,alpha=0.7,palette=custom_palette1,ax=ax[1],zorder=0)
ax[1].legend(bbox_to_anchor=(2,1),loc='center right')
ax[1].set_xlim(0,1.2) #lim for flow velocity
#### end of older plotting of sensitivity #####



#### here lab 85 to see sensitivity as settling to breakup ratio fig.3###
path85='/gpfs/work/lie/run/par/4lab85/range/bale2002.nc' ##RMSES,old,range
path780='/gpfs/work/lie/run/par/4lab780/0.1RMSE1range/bale2002.nc' ##RMSES/' 


##fettref='/gpfs/work/lie/run/par/6oszmay/testplot/bale2002_.nc'

def ratiocalculation(path,labelinput,c):
	ncvi=netCDF4.Dataset(path).variables
	ratiocali=-squeeze(ncvi['agg_sinkinglpm'][fr:to,lfr:lto])/squeeze(ncvi['agg_Breakup'][fr:to,lfr:lto]) #fr:to
	G4h = squeeze(ncvi['agg_G'][fr:to,lfr:lto])#.astype(bool)  #
	ax.scatter(G4h,ratiocali,alpha=0.5,s=50,label=labelinput,color=c) #label=labelinput,
	return ratiocali



plt.close('all')
print ('plot settling to breakup ratio of lab85 test')
#################### fig3 the 3 cases settling to breakup ratio under different G during one simulated tidal cycle###
plt.close('all')
fig, ax = plt.subplots()
ax.axhline(y=1, color=".7", dashes=(2, 1),lw=0.5, zorder=0)
#ratiocalculation(path780)
ratiocalculation(baleref,'high SPMC, h=0.28m (lab)','k') #'lab3500'
ratiocalculation(path85,'low SPMC, h=0.28m (lab)','grey') #'lab85'
ratiocalculation(fettref,'low SPMC, h=20m (field)','dodgerblue') #'field20'
x=np.logspace(-3,1.1,100)
y=x**-1.5
plt.plot(x,y,color='k',ls='--', lw=1,label='slope=-1.5')

plt.yscale('log')
plt.xscale('log')
ax.set_xlim(0,None)
ax.set_xlabel('turbulence $G$ (s$^{-1}$)')
ax.set_ylabel('settling-to-breakup ratio')
#plt.legend()

ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.2),ncol=2)


###################### end of fig3 3 slopes of settling to breakup ratio ##############
plt.show()
exit()


ncvi=netCDF4.Dataset(path85).variables
ratiocali=-squeeze(ncvi['agg_sinkinglpm'][fr:to,lfr:lto])/squeeze(ncvi['agg_Breakup'][fr:to,lfr:lto]) #fr:to
G4h = squeeze(ncvi['agg_G'][fr:to,lfr:lto])#.astype(bool)  #
plt.close('all')
print ('plot settling to breakup ratio of lab85 test')
fig, ax = plt.subplots()
ax.scatter(G4h,ratiocali)
plt.yscale('log')


plt.close('all')
g=sns.relplot(
    data=df, x=xplot,y=yplot,kind='line',hue='Amplitude ref',lw=5,alpha=0.7,palette=custom_palette #mako ,hue="depth(m)",style="Amp",palette=sns.color_palette("viridis_r", as_cmap=True)mako  palette=sns.color_palette('rocket',n_colors=len(filelist)
)#lw=5,
(g.map(plt.axhline, y=1, color=".7", dashes=(2, 1),lw=0.5, zorder=0).set_axis_labels(xplot,yplot).tight_layout(w_pad=0))

ax=g.axes[0,0]
#ax2=ax.twiny() #not successful to get double xlabel also for turbulence as high non-linearity
##ticklist=df['G (s$^{-1}$)'].tolist()
#ax2.set_xticks([0,0.4,0.8,1.2]) #set ticks for non-linear scale for G as secondary axis.
#ax2.plot(df['velocity'],df['ratio']) #df['G (s$^{-1}$)']
##[0.1,0.5,0.9,1.3] for lab
##[0.01,0.02,0.03,0.04]
#ax.set_xlim(0,1.2) #lim for fow velocity
plt.yscale('log')
g.set(ylim=(0,None))
#plt.xscale('log')
plt.show()
exit()








################# end of plotting ###################

