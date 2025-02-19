#!/usr/bin/python3
####code for evaluation of the slurm job, only dealing with txt files#############
#############for field case of 6 parameters###########################################
##tmzjulnew0

import pandas as pd
import matplotlib
from matplotlib import rc
from matplotlib import ticker, cm
import netCDF4
from pylab import *
import sys
#import cftime #netcdftime
import numpy
import pickle
from mpl_toolkits import mplot3d
import seaborn as sns
from collections import defaultdict
import scipy
from heapq import nsmallest
import subprocess


####check if drop the runs that was contant!!!!!!!!!!!!!############

#from astropy.visualization import simple_norm
#print('matplotlib: {}'.format(matplotlib.__version__))

pvar0 ='v1' #input("Enter xaxis for plotting (v?): ") #'v0'

####log file for parameters####### tmzjulnew0
dfp = pd.read_csv('/gpfs/work/lie/results/oszmaynew3/current.log', sep="-",header=None)#, names=['id','param','value']) #parameters       run/test/ bigrun numpy.linalg.LinAlgError: 2-th leading minor of the array is not positive definite
#currentsample.log
dfp.columns = ['id','param','value']
sdfp = dfp.sort_values(by=['id'], ascending=True) #sorted
#idp0 = sdfp['id'].values[0]
idp=sdfp['id'].to_numpy()
p0raw=sdfp.loc[sdfp['param'] == 'breakup_factor']  #extract all parameter 0
p1raw=sdfp.loc[sdfp['param'] == 'coagulation_rate'] #
p2raw=sdfp.loc[sdfp['param'] == 'ks']   #'dens_lpm'
###add more parameter####
p3raw=sdfp.loc[sdfp['param'] == 'fractal_dimension1'] 
p4raw=sdfp.loc[sdfp['param'] == 'pws'] # 
p5raw=sdfp.loc[sdfp['param'] == 'kws'] # 
p6raw=sdfp.loc[sdfp['param'] == 'rho'] #'kd'

###checking raw values of parameters#######
rawlist=[p0raw,p1raw,p2raw,p3raw,p4raw,p5raw,p6raw]
print ("p6 values are:",p6raw['value'])


###check raw list if the parameters are the same value########
#namelist=['breakup_factor','coagulation_rate','ks','fractal_dimension1','kbottom','tauc_const','Xsize'] # 
namelist=['breakup_factor','coagulation_rate','ks','fractal_dimension1','pws','kws','rho'] #'Xsize' , 'kbottom','tauc_const',
print ('first checking if variations are done')
for i in range(0,7):
	xi=sdfp.loc[sdfp['param'] == namelist[i]] 
	yi=xi['value']
	#yi.to_string(index=False)
	yi.values.tolist() 
	#print (len(yi)) #yi[0],yi[10]
	tof=len(set(yi)) <= 1
	if tof==True:
		print ("param nr.",i,namelist[i], ":all same value:",tof,"!!!!!")  #yi.count(yi[0]) == len(yi),
	else:
		print ("param nr.",i,namelist[i], "worked")


####relative error result file#############
dfe = pd.read_csv('/gpfs/work/lie/results/oszmaynew3/relativeerror.txt', sep=" ",header=None)#, names=['id','v2','v3','v4','v5','v6','v7']) #errors bigrun
#relativesample.txt

#dfe.columns = ['id','mse','mses','msesum','valley','peaknr','errsum']
#dfe.columns = ['id','mse','mses','msesum','valley','peaknr','errsum']

#dfe.columns = ['id','mse','mses','msesum','valley','peaknr','err_dmax','err_dmin','err_dmean','err_drange','err_cmax','err_cmin','err_cmean','err_crange','rangesum','errsum']
#dfe.columns = ['id','mse','mses','msesum','valley','peaknr','err_dmax','err_dmin','err_dmean','err_drange','err_cmax','err_cmin','err_cmean','err_crange','rangesum','stdd','stdc','errsum']
dfe.columns = ['id','mse','mses','msebt','msesdt','peaknr','err_dmax','err_dmin','err_dmean','err_drange','err_cmax','err_cmin','err_cmean','err_crange','rangesum','stdd','stdc','errsum']


sdfe = dfe.sort_values(by=['id'], ascending=True) #sorted
ide=sdfe['id'].to_numpy()
ide=ide.astype(int)
print ("valid runs are:",ide)

####add for plotting##
#sdfe['errsum_r']=1+1/sdfe['errsum'] #just for plotting  #abs(1-sdfe['valley'])+sdfe['msesum']
#plt.plot(sdfe['errsum'],sdfe['errsum_r']) #1+1/sdfe['errsum']
#plt.show()
#exit()
######drop duplicated and get the best runs #####
####filtering of run results######!!!!!!!!!
simple=sdfe.drop_duplicates(subset=['errsum'], keep='last')
#######recalculating range error as err_range, nomalized by mean value#######################
simple['err_cr']=(simple['err_cmax']+simple['err_cmin'])#/40
simple['err_dr']=(simple['err_dmax']+simple['err_dmin'])#/80 #200
simple['err_range']= simple['err_cr']+simple['err_dr']
simple['msesum']=simple['mse']+simple['mses']+simple['err_range']

simple['sum']=simple['mses'] #simple['mse']/40+simple['mses']/80 #simple['err_cr']*0.1+simple['err_dr'] #
#+0.5*(simple['err_range'])
#print (simple)
#exit()

#bestinput=5000 #10000 #40000 #20
bestinput=round(0.001*len(simple))  #simple1 for lab
print ("best 1% is _ many runs", bestinput)

#########here for peak-to-valley for lab one #################################
#bestv=simple[(simple['valley']<=1) & (simple['valley']>=0)] #filter reasonable valley range (should be [0,1])
#bestn=bestv[(bestv['peaknr']==8)] #==8  filter peaknr=8 for 2 peak shape
#best=bestn.nsmallest(bestinput,'errsum') #errsum_r .index.get_level_values(1) #according to error sum best runs 
##worse=bestn.nlargest(bestinput,'errsum') #errsum_r .index.get_level_values(1) #according to error sum best runs 


####here for field case##########################3
simple1=simple[(simple['peaknr']>2)] #filter peaknr>1 for non constant numbers
#best=bestv.nlargest(10,'valley')
#bestmse=bestv.nsmallest(100,'msesum')   
#bestmse=simple1.nsmallest(bestinput,'msesum')   

bestspmc=simple1.nsmallest(bestinput,'mse')   
bestsize=simple1.nsmallest(bestinput,'mses')   

bestspmc_t=simple1.nsmallest(bestinput,'msebt')   
bestsize_t=simple1.nsmallest(bestinput,'msesdt')   

bestrange=simple1.nsmallest(bestinput,'rangesum') 

beststdc=simple1['stdc'].sort_values(ascending=True).head(bestinput)
#simple1.nsmallest(bestinput,'stdc') 

beststdd=simple1['stdd'].sort_values(ascending=True).head(bestinput)
#simple1.nsmallest(bestinput,'stdd') 

bestmsesum=simple1.nsmallest(bestinput,'msesum')   

bestsum=simple1.nsmallest(bestinput,'sum') #trying matrix!!!!!!!!!!!!!!!!!!!!!!1
 
#therun=bestn.nsmallest(1,'errsum' ) #select the best 1 run already filtered peaks, for lab

#therun=simple.nsmallest(1,'msesum')

#print ("THE run is number: ")
#print (theid)

####now try intercept################
#inter=bestsum #bestspmc # bestsize#simple1 #simple
inter=bestsum[['id','sum']] ###only keep id and sum to make smaller dataframe!!!!!!!!!!!!!

#inter=bestspmc[bestspmc['id']&bestsize['id']]
#inter=simple1[bestspmc['id']&bestsize['id']] #bestspmc['id']&bestsize['id'] #altogether with range??? 
#inter=bestspmc[bestspmc['id']&bestsize['id']&bestrange['id']&beststdc['id']&beststdd['id']] ##altogether with range??? 
#inter=bestspmc[bestspmc['id']&bestsize['id']]  #bestspmc #  bestsize['id']    best[best['id']&bestmse['id']] #overlapping between two best schemes
print (inter) ###here!
idinter=inter['id']  #"sorting according to best 10 total & intersection of best 100 MSE"

############define which error to plot#######################################################
ploterr='sum' #'mses' # 'msesum' #'err_range' #'err_drange' #'msesdt' #'stdd' #'msesum'  #'rangesum' #'mses' #'stdd'
print ("error plotted here is", ploterr)
therun=inter.nsmallest(1,ploterr) #ploterr #defining the best run
theid=therun['id']
print (therun)
#exit()
#print ("best "+str(bestinput)+" runs used")  #all runs recorded are used


#test=pd.DataFrame()
test4=pd.DataFrame()
p0=pd.DataFrame()
p1=pd.DataFrame()
p2=pd.DataFrame()
p3=pd.DataFrame()
p4=pd.DataFrame()
p5=pd.DataFrame()
p6=pd.DataFrame()

#print(p1raw)
#exit()

def dffilter(x):
	p0x=pd.DataFrame()
	p1x=pd.DataFrame()
	p2x=pd.DataFrame()
	p3x=pd.DataFrame()
	p4x=pd.DataFrame()
	p5x=pd.DataFrame()
	p6x=pd.DataFrame()
	for i in x[:]: # ide[:]:##idinter: #idbest[:]:#  ##ide[:]: #loop to get parameter
   		locp0x=p0raw.loc[p0raw['id']==i]
   		p0x=pd.concat([p0x,locp0x])
   		locp1x=p1raw.loc[p1raw['id']==i]
   		p1x=pd.concat([p1x,locp1x])
   		locp2x=p2raw.loc[p2raw['id']==i]
   		p2x=pd.concat([p2x,locp2x])
   		locp3x=p3raw.loc[p3raw['id']==i]
   		p3x=pd.concat([p3x,locp3x])
   		locp4x=p4raw.loc[p4raw['id']==i]
   		p4x=pd.concat([p4x,locp4x])
   		locp5x=p5raw.loc[p5raw['id']==i]
   		p5x=pd.concat([p5x,locp5x])
   		locp6x=p6raw.loc[p6raw['id']==i]
   		p6x=pd.concat([p6x,locp6x])
	test1x=p0x.merge(p1x,on='id',how='inner')
	testx=test1x.merge(p2x,on='id',how='inner')
	test2x=testx.merge(p3x,on='id',how='inner')
	test3x=test2x.merge(p4x,on='id',how='inner')
	test4x=test3x.merge(p5x,on='id',how='inner')
	test4x.columns=['id','p0','v0','p1','v1','p2','v2','p3','v3','p4','v4','p5','v5']  #rename	#test5x=test4x.merge(p6x,on='id',how='inner')
#	test5x.columns=['id','p0','v0','p1','v1','p2','v2','p3','v3','p4','v4','p5','v5','p6','v6']  #rename
#	test2x.columns=['id','p0','v0','p1','v1','p2','v2','p3','v3']
	datax=test4x.merge(simple1,on='id',how='inner') #test   #full frame  #test for 3 params	
#	datax=test2x.merge(simple,on='id',how='inner') #test   #full frame  #test for 3 params	  #sdfe is simple without dropping duplicates
	return datax


#pvar1 = input("yaxis (v?): ")
#print("You entered: " + pvar0, pvar1)

names=['v0','v1','v2','v3','v4','v5','v6']  ##############!!!!!

#nruns=1000
#ids=bestn.nsmallest(nruns,'errsum')
###get the runs with best nruns################3
#datarf=dffilter(ids['id']) #datarefiltered only for valley
#print ("dataframe used best "+str(nruns)+" runs")

######datarf all the calculation later are done################################
datarf=dffilter(bestsum['id'][:]) #simple1 idinter[:]  simple1['id'][:]  bestmse['id'][:]  (idbest[:) #best for smallest errsum  #idinter[:] ##datarf is filtered with valley and best n values
print ("dataframe created for plotting")
print (datarf)

#exit()


translation = {39: None}

####find the median for comparision#########################
medians=[] 
medians.append(000)
paramedian=[]
for i,n in enumerate(names[0:6]): #[0:7] #namelist is the long original name, names are v0, v1...
	#x=namelist[i]
	#print (eval(x))
	xmax=datarf[n].max()
	xmin=datarf[n].min()
	xmiddle=(xmax+xmin)/2
	medians.append(namelist[i])
	medians.append(xmiddle) #datarf[n].median()
	paramedian.append(xmiddle)#datarf[n].median()
print ("median value of each column")
clear=str(medians).translate(translation)  #remove apostroph/ single quotes 

print ("median is",medians)
print ("median string",clear)
f=open('/gpfs/work/lie/results/medianoszmay.txt','w') ##'ab' 'w'
for i in medians:
  f.write(str(i)+'  ')#'\n'
#f.write(b'\n')
f.close()

#print ("paramedians are",paramedian)
##get the parameter for the best run#########3
for j in theid[:]:
	#therunpara=np.where(data['id']==j)
	therunpara=datarf[datarf['id'] == j]
print (therunpara, "best run id is: ", j)
np.savetxt('therunoszmay.txt', therunpara, fmt='%s')
np.savetxt('paramtochangeoszmay.txt', therunpara, fmt='%s')
#print (list(data)) 

########added plotting the time series, call another python script############################################
print ("start plotting time series with data")
#subprocess.call("/gpfs/work/lie/bale2002/parameterauto6.sh",shell=True)  #!!!a different one!
print ("end of parameterauto.sh")

##############################################################################################################
def findparameter(x):
	runnr=datarf[datarf['id'] == x]
#	print (runnr, "best run id is: ",x)
	f=open('/gpfs/work/lie/results/paramtochangeoszmay.txt','ab')
	np.savetxt(f,runnr,fmt='%s',delimiter=' ')
	f.close()
	return runnr


def finderror(x,y):
	#runnr=datarf[datarf['errsum'] == x]
	runnr=datarf[datarf[y] == x]
#	print (runnr, "best run id is: ",x)
	f=open('/gpfs/work/lie/results/paramtochangeoszmay.txt','ab')
	np.savetxt(f,runnr,fmt='%s',delimiter=' ')
	f.write(b'\n')
	f.close()
#	np.savetxt('paramtochange6.txt', runnr, fmt='%s')  #runnr'+str(x)+
	return runnr


#print (findparameter(15966))
#exit()

##############averaging the parameters for best 1% runs ##########################
meanpara=[]
stdpara=[]
kb=[]
kc=[]
for i,n in enumerate(names[0:6]): #[0:7] #namelist is the long original name, names are v0, v1...
	#x=namelist[i]
	#print (eval(x))
	#meanpara.append(namelist[i])
	#meanpara.append(datarf[n].median())
	meanpara.append(datarf[n].mean())
	stdpara.append(datarf[n].std())
	if i==0:
		kb.append(datarf[n])
	if i==1:
		kc.append(datarf[n])
print ("mean value of each column")
print (meanpara)
print ("std value of each column")
print (stdpara)

rkb = [a_elem // b_elem for a_elem, b_elem in zip(kb,kc)]
print ("kb/kc")
print (mean(rkb),std(rkb))
exit()
####here!!!!###
####find certain runs, least error at certain parameter combo####
#rungroup=datarf.groupby(['v3','v1','id']).agg({ploterr:'min'}).reset_index()
run00=datarf[datarf['v5'] == 2700]

#run10=run00[run00['v1']==90].nsmallest(1,ploterr)#.astype(int)
#run30=run00[run00['v1']==100].nsmallest(1,ploterr)#.astype(int)
run110=run00[run00['v1']==110].nsmallest(1,ploterr)#.astype(int) #.to_numpy()
run120=run00[run00['v1']==120].nsmallest(1,ploterr)#.astype(int) #.to_numpy()
run130=run00[run00['v1']==130].nsmallest(1,ploterr)#.astype(int) #.to_numpy()

print ("id for runs with v1 = 110")
print (run110['id'],run110[ploterr])
print ("id for runs with v1 = 120")
print (run120['id'],run120[ploterr])
print ("id for runs with v1 = 130")
print (run130['id'],run130[ploterr])
#print ("id for runs with v1 = 70")
#print (run70['id'],run70[ploterr])
#print ("id for runs with v1 = 90")
#print (run90['id'],run90[ploterr])
#exit()

n00=19684 #14818
n11=19245 #10843
n22=19637 #8880 #8931 #8930
findparameter(n00)
findparameter(n11) #10842
findparameter(n22) #8880

#exit()
kc110=datarf[datarf['id'] == n00]
kc120=datarf[datarf['id'] == n11] #10842
kc130=datarf[datarf['id'] == n22]

#exit()

#####2D plot xaxis one parameter, y another, with scatter of min error##############
matplotlib.rcParams.update({'font.size': 22})
fig,axs = plt.subplots(2,3,constrained_layout=True,sharex=True) #2,3,
axs = axs.ravel()

#xt2min = datarf.groupby([pvar0])['v1','errsum'].min().reset_index()
#print (xt2min)
#plt.scatter(xt2min[pvar0],xt2min['v1'],c=xt2min['errsum'])
#plt.show()

###################errsum for lab, others testing#############


for i,n in enumerate(names[0:6]): #[2:6] [0:7]
	print (n)
	#datac=data[[pvar0,t,'errsum']]
	if n == pvar0:
		axs[i].plot(datarf[pvar0],datarf[pvar0])
		continue
	else:
		#xt2min = datarf.groupby([pvar0])[n,'errsum'].min().reset_index()
		xt2min = datarf.groupby([pvar0, n]).agg({ploterr:'min'}).reset_index()
		#print (xt2min)
		im=axs[i].tricontourf(xt2min[pvar0],xt2min[n],xt2min[ploterr],levels=50,cmap='YlGn_r', extend='max')#/xt2min[ploterr].max() ,, cmap='viridis_r')
#		im=axs[i].scatter(xt2min[pvar0],xt2min[n],c=xt2min[ploterr],cmap='viridis_r')  #scatter plot

#		xre = scipy.ndimage.zoom(xtest[pvar0], 3) #even more refined!!!
#		yre = scipy.ndimage.zoom(xtest[t], 3)
#		zre = scipy.ndimage.zoom(xtest['errsum'], 3)
#		im=axs[i].tricontourf(xre,yre,zre,levels=14,cmap='viridis_r') #interpolate data
		print (xt2min[n],xt2min[ploterr])
		axs[i].scatter(therunpara[pvar0],therunpara[n],marker='^', c='red', s=30)
		#axs[i].set(xlim=(0,5000))
		axs[i].set_title(n,fontsize='small',  loc='left') #
		axs[i].scatter(kc110[pvar0],kc110[n],marker='o',c='orange',s=50) #c='blue',
		axs[i].scatter(kc120[pvar0],kc120[n],marker='+',c='k',s=50) #c='blue',
		axs[i].scatter(kc130[pvar0],kc130[n],marker='.',c='blue',s=50) #c='blue',
#axs[3].set_facecolor('black')  ##set back ground color 
###scatter of median values########
#axs[0].scatter(paramedian[3],paramedian[2],marker='o', c='orange', s=30)
#axs[2].scatter(paramedian[3],paramedian[4],marker='o', c='orange', s=30)
#axs[3].scatter(paramedian[3],paramedian[5],marker='o', c='orange', s=30)
fig.supxlabel(pvar0)
cbar=fig.colorbar(im,ax=axs.flat)
cbar.ax.get_yaxis().labelpad = 18
cbar.ax.set_ylabel(ploterr, rotation=270)

plt.show()
exit()

#a=findparameter(14357)
#findparameter(40661)
#print (a['errsum'])
#findparameter(35126)




#for i,n in enumerate(np.unique(datarf[pvar0])[:]):  #names[:]
#	print (n)

	
#	datar=datarf[['id',t,'errsum']]
	#vu=np.unique(datac[pvar0])
	#for ii in vu:
#	x2tmin = datar.groupby([t]).min().reset_index()
#	x2tmean=datar.groupby([t]).mean().reset_index()
#	axs[i].scatter(x2tmin[t],x2tmin['errsum'],label='min error')
#	axs[i].scatter(x2tmean[t],x2tmean['errsum'],label='mean error best'+str(nruns))
#	axs[i].set_title(t,fontsize='small',  loc='left') #
#lines, labels = fig.axes[-1].get_legend_handles_labels()   
#fig.legend(lines, labels, loc='upper right',) #'upper right'
#fig.supxlabel(pvar0)
#plt.show()

##############################################3
#print(datarf[datarf['id'] == 12220])
#datarf is correctly constructed!
#exit()


########1D plot for each value of certain parameter###################
######1 line for min error, 1 line for min error (0.1%)###############
matplotlib.rcParams.update({'font.size': 22})

close('all')
#matplotlib.use('Qt5Agg')
fig,axs = plt.subplots(2,4) #,constrained_layout=True,sharex=True
axs = axs.ravel()
howmany=2 #how many least errorsum runs in each parameter
#xtest=datarf
#print (", " .join(names))
#xtest=datarf.groupby(names[:])[['errsum','msesum','id','v0', 'v1', 'v2', 'v3', 'v4', 'v5']].agg('min')  #", " .join(names)
#xtest = datarf.groupby(names[:])[['errsum','msesum','id','v0', 'v1', 'v2', 'v3', 'v4', 'v5']].min()
#xtminreal=xtest.nsmallest(howmany,'errsum') #wrong only 2 
####find smallest error for all parameter values############
besthowmany=2 #3
for i,j in enumerate(names[:]): #for each parameter
	para=[]
	err=[]
	n2find=[]
	print ('for',j)
	for ii in np.unique(datarf[j]): #for each value of this parameter
		thisdf = datarf[datarf[j] == ii] #
		MIN=thisdf[[ploterr,'id']].min().min()  #,j with j mess up with err #'errsum'
		theline=finderror(MIN,ploterr)       ###########???????here changed to ploterr to stay consistent 
		para.append(ii)
		err.append(MIN)
		n2find.append(theline['id'].tolist()) #.tolist() 
		#print (n2find) #
		#print(ii,MIN)
	test=nsmallest(besthowmany, err)#[-1] ##smallest error run for certain parameter
	test1=[err.index(testi) for testi in test]  #find its relative place
	print(test, test1)
	#testtest=''.join(map(str,test1)) #[1:-1]) change to suitable list
#	for testtesti in testtest:
#		print (testtesti)
#		z=n2find[int(testtesti)]
#		print (z)
	n2test=[ n2find[test1i] for test1i in test1 ] #find its ID  #[int(test1i)] 
	#n2test=''.join(map(str,n2test))
	for jj in n2test:  
		for jjj in jj:
			findparameter(jjj) #get the parameter to file for best n runs
			print (jjj)
	#n2find=np.array(n2find)
	#y=[ findparameter(item) for item in n2test ] #find for each value each parameter..
	
	axs[i].scatter(np.array(para),np.array(err))#,label='mean error best'+str(bestinput)) #nruns
plt.show()

exit()


#axs[i].set_title(t,fontsize='small',  loc='left') #
#lines, labels = fig.axes[-1].get_legend_handles_labels()   
#fig.legend(lines, labels, loc='upper right',) #'upper right'

####################################################################end##############################################
####older two set plot scatter######correct plot wrong run numbers######
fig,axs = plt.subplots(2,3) #,constrained_layout=True,sharex=True
axs = axs.ravel()

for i,t in enumerate(names[:]):
	print (t)
#	datar=datarf[['id',t,'mse','mses','errsum']]
	datar=datarf
	#vu=np.unique(datac[pvar0])
	#for ii in vu:
#	xtmin = datar.groupby(t)['errsum'].min().reset_index()

####both solution works###################
######correct plotting but incorrect number of run recording#############
######always recording the smallest errors###########################

	#xtmin = datar.groupby(t)[['errsum','id',t]].min() #.min()#.reset_index()  (t,as_index=False)
	xtmin = datar.groupby(t)[['errsum','rangesum','id',t]].agg('min')   #'msesum'

	print (xtmin)
	#xtminreal=xtmin.groupby('errsum').min()
	xtminreal=xtmin.nsmallest(howmany,'errsum')


#	xtmin1 = datar.groupby(t).agg({'errsum': ['min','max']})

	#xtmin1 = datar.groupby(t).agg({'errsum': ['min']})  #'mean', 'min', 'max'
	xtmean=datar.groupby(t)['errsum'].mean().reset_index()
	axs[i].scatter(xtmin[t],xtmin['errsum'],label='min error')
#	print (xtmin[['id',t,'errsum']]) # here the min error is not correct!!!!!!!

	print (xtminreal) # here the min error is not correct with run numbers !!!!!!!
#	print (xtmin1['errsum'])


#	testindex=xtmin.loc[xtmin.groupby('errsum')['id'].idxmin()].reset_index(drop=True) #sort lowest 2 values

	#ntofind=[]
	#ntofind.append(testindex.loc[0:1]['id'])

#	n2find=xtminreal['id'].tolist() 
#	print (n2find) #

#	n2find=testindex.loc[0:1]['id'].tolist() #print the smallest error and correponded run number

#	y=[ findparameter(j) for j in n2find[:] ] #findparameter(j)  print(j) 

	#print (xtmin1)
	axs[i].scatter(xtmean[t],xtmean['errsum'],label='mean error best'+str(bestinput)) #nruns
	axs[i].set_title(t,fontsize='small',  loc='left') #
lines, labels = fig.axes[-1].get_legend_handles_labels()   
fig.legend(lines, labels, loc='upper right',) #'upper right'
#fig.supxlabel(pvar0)

plt.show()

#tester=datarf[datarf['id'] == 14357] #################try to fix 
#datarf[datarf['id'] == 40661]
#print(tester) #[['id','errsum']]



exit()










data=dffilter(idbest)


###define the best or all the results############ide is for all run###

#for i in idbest[:]: # ide[:]:##idinter: #idbest[:]:#  ##ide[:]: #loop to get parameter
#   locp0=p0raw.loc[p0raw['id']==i]
#   p0=pd.concat([p0,locp0]) 

#   locp1=p1raw.loc[p1raw['id']==i]
#   p1=pd.concat([p1,locp1]) 

#   locp2=p2raw.loc[p2raw['id']==i]
#   p2=pd.concat([p2,locp2]) 

#   locp3=p3raw.loc[p3raw['id']==i]
#   p3=pd.concat([p3,locp3]) 

#   locp4=p4raw.loc[p4raw['id']==i]
#   p4=pd.concat([p4,locp4]) 

#   locp5=p5raw.loc[p5raw['id']==i]
#   p5=pd.concat([p5,locp5]) 

#test1=p0.merge(p1,on='id',how='inner')
#test=test1.merge(p2,on='id',how='inner')
##add more params##
#test2=test.merge(p3,on='id',how='inner')
#test3=test2.merge(p4,on='id',how='inner')
#test4=test3.merge(p5,on='id',how='inner')
#test4.columns=['id','p0','v0','p1','v1','p2','v2','p3','v3','p4','v4','p5','v5']  #rename
##test.columns=['id','p0','v0','p1','v1','p2','v2']  #rename
#####
#data=test4.merge(sdfe,on='id',how='inner') #test   #full frame  #test for 3 params
## 'id','p0','v0','p1','v1','p2','v2','p3','v3','p4','v4','p5','v5','mse','mses','msesum','valley','peaknr','errsum','errsum_r'
##  id param_x  value_x  param_y  value_y param  value   mse  mses msesum    valley  peaknr  errsum   errsum_r


#exit()
#print (data)
#print (data['v3'])

#print (data['value_y'])

###plotting one parameter to error#######
plt.rcParams.update({'font.size': 22})
fig, ax = plt.subplots()
ax.set_xlabel('coagulation') #'breakup'#'coagulation'
ax.set_ylabel('error sum')

#ax.scatter(p0,sdfe['errsum'].iloc[:],label='breakup') #plotting id and first parameter

#ax.scatter(p0['value'],abs(1-sdfe['valley'].iloc[:]),color='blue',label='breakup') #plotting id and first parameter
#ax.scatter(p0['value'],best['errsum_r'] .iloc[:],color='blue',label='breakup') 

#ax.scatter(p1['value',sdfe['errsum'].iloc[:],color='orange',label='coagulation') 
#ax.scatter(p1['value'],abs(1-sdfe['valley'].iloc[:])+sdfe['msesum'].iloc[:],color='orange',label='coagulation')

#ax.scatter(p2['value'],sdfe['errsum'].iloc[:],color='green',label='sinking') 
#ax.scatter(p2['value'],abs(1-sdfe['valley'].iloc[:])+sdfe['msesum'].iloc[:],color='green',label='sinking') 

###1-valley: abs(1-sdfe['valley'].iloc[:])
###total error: abs(1-sdfe['valley'].iloc[:])+sdfe['msesum'].iloc[:]

ax.scatter(data['v1'],data['mse'],color='blue',label='v1')  #value_x  breakup
plt.legend()
#print(len(sdfe['errsum'].iloc[:]))
print ("length of parameter 0 is",len(p0))

#exit()

x=data['v2']         #'v3'  ['value'] #'value_x'
##compound plot of 4 error matrix###
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=False) #,constrained_layout=True
ax1.scatter(x,data['mse'],c='blue')#,label='MSE_SPMC')
ax1.set_title('MSE_SPMC',loc='left')
ax2.scatter(x,data['mses'],c='orange')#,label='MSE_D50')
ax2.set_title('MSE_D50',loc='left')
ax3.scatter(x,abs(1-data['valley']),c='green')#,label='Valley')
ax3.set_title('Valley',loc='left')
ax4.scatter(x,data['errsum_r'],c='k')#,label='Total error')
ax4.set_title('Total error',loc='left')
fig.text(0.5, 0.05,'Sinking parameter' , ha='center')  #Sinking parameter  "fractal dimension"
plt.tight_layout()
#plt.legend()

###scatter plot for parameter and error #####
fig,ax = plt.subplots()
#ax = plt.axes(projection='3d')
#E1=data['value'] #p2['value']#sdfe['errsum_r'] #sdfe best log(abs(1-sdfe['valley'].iloc[:])+sdfe['msesum'].iloc[:])
p=ax.scatter(data['v0'],data['v1'],c=data['errsum'], s=30, cmap='viridis') #'BuGn_r' best marked with dark green
#p=ax.scatter(p0['value'], p1['value'], c=E1, s=30, cmap='viridis') #'BuGn_r' best marked with dark green

#x,y=np.meshgrid(p0['value'], p1['value'])
#z = np.array(sdfe['errsum_r']) #p2['value']
#z = z.reshape((len(x), len(y)))
#plt.contourf(x,y,z)

ax.set_xlabel('breakup') 
ax.set_ylabel('coagulation') #''
#ax.set_zlabel('errorsum')
#plt.colorbar(p)
cbar = plt.colorbar(p)
cbar.ax.set_ylabel('error_sum', labelpad=20,rotation=270)

#print (E1)

#plt.scatter(sdfe['id'],sdfe['errsum'])



#####contour plot has multiple entry, plot the min ###############

fig,axs = plt.subplots(2,3,constrained_layout=True) #,sharex=True
axs = axs.ravel()
#datap=data[['v0','v1','errsum']]
#x1=np.unique(datap['v0'])

#xtest = datap.groupby(['v0','v1']).mean().reset_index() #
#z=xtest['errsum']
#ax.scatter(xtest['v0'],xtest['v1'],c=z,cmap='Reds')
#plt.tricontourf(xtest['v0'],xtest['v1'],z,levels=14,cmps='Reds')
#plt.show()





#######tricontour plot for all best runs#################
#######not enough then scatter###################
for i,t in enumerate(names[:]):
	print (t)
	datac=data[[pvar0,t,'errsum']]
	if t== pvar0:
		axs[i].plot(datac[pvar0],datac[pvar0])
		continue
	else:
		dataci=datac['errsum'].values.tolist() 
		datact=datac[t].values.tolist() 
		datacp=datac[pvar0].values.tolist() 
		#print (len(yi)) #yi[0],yi[10]
		if len(set(dataci)) <= 1 or len(set(datact)) <= 1 or len(set(datacp)) <= 1: #add loop to scatter when not enough point to contour
			print (t,"all values are same")
			im=axs[i].scatter(datac[pvar0],datac[t],c=datac['errsum'],cmap='viridis_r')
			print ()
		else:
			xtest = datac.groupby([pvar0,t]).min().reset_index() #mean()  make average of error 
			im=axs[i].tricontourf(xtest[pvar0],xtest[t],xtest['errsum'],levels=14,cmap='viridis_r')
#		xre = scipy.ndimage.zoom(xtest[pvar0], 3) #even more refined!!!
#		yre = scipy.ndimage.zoom(xtest[t], 3)
#		zre = scipy.ndimage.zoom(xtest['errsum'], 3)
#		im=axs[i].tricontourf(xre,yre,zre,levels=14,cmap='viridis_r') #interpolate data

		axs[i].scatter(therunpara[pvar0],therunpara[t],marker='^', c='red', s=30)
		#axs[i].set(xlim=(0,5000))
		axs[i].set_title(t,fontsize='small',  loc='left') #
fig.supxlabel(pvar0)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#fig.colorbar(im, cax=cbar_ax)
cbar=fig.colorbar(im, ax=axs.flat)
cbar.ax.get_yaxis().labelpad = 18
cbar.ax.set_ylabel('errsum', rotation=270)
plt.show()
exit()
#print (therunpara[pvar0])

def densityplot(xin,yin,zin,**plt_kwargs):
	#datap=data[[xin,yin,zin]]
	sns.displot(data,x=xin,y=yin, kind='kde',cmap='mako_r',fill=True) #,cut=0,

#	norm = plt.Normalize(data[zin].max(),data[zin].min())
#	sm = plt.cm.ScalarMappable(cmap='mako_r', norm=norm) #
#	sm.set_array([])
#	plt.colorbar(sm)

#	sms=plt.scatter(data[xin],data[yin],c=data[zin], s=30, cmap='Reds') #'BuGn_r' best marked with dark green
#	plt.colorbar(sms)
	plt.scatter(therunpara[xin],therunpara[yin],c='red', s=30)
	return (xin, yin, zin)

list_im=[]
for i,t in enumerate(names[:]):
	print (t)
	if t==pvar0:
		continue
	else:
		densityplot(pvar0,t,'errsum')
		plt.savefig(str(pvar0) + '-' + str(t) +'.png')
		list_im.append(str(pvar0) + '-' + str(t) +'.png')

######combining plot##########
#listimage=["v0-v1.png","v0-v2.png","v0-v3.png","v0-v4.png"]
import PIL
from PIL import Image

#list_im = listimage
imgs    = [ PIL.Image.open(i) for i in list_im ]
# pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )

# save that beautiful picture
imgs_comb = PIL.Image.fromarray( imgs_comb)
imgs_comb.save( 'Trifecta.png' )    
img = Image.open('Trifecta.png')
img.show() 

exit()

## for a vertical stacking it is simple: use vstack
#imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
#imgs_comb = PIL.Image.fromarray( imgs_comb)
#imgs_comb.save( 'Trifecta_vertical.png' )




###trying to subplot from loop## cannot reshape z#########
#for t, ax in zip(names, axes.flatten()): ##########here!!!!!!!!!!################??????????
#    sns.boxplot(y=name, x= "a", data=df, orient='v', ax=ax)
#for i,t in enumerate(names):
#    	sns.boxplot(y=t, x= "a", data=df, orient='v', ax=axs[i % 2], palette=flatui)
#	datap=data[[pvar0, t,'errsum']] #  'msesum'
#	datap=pd.DataFrame(datap,index=[0])
#	print (datap)
#	print ("pvar0 is:",datap[pvar0], "t is:",datap[t])
#	sns.displot(datap,x=pvar0, y=t,kind='kde',orient='v')
#	sns.boxplot(y=t, x= pvar0, data=datap, orient='v', ax=ax)
#	sns.displot(datap,x=pvar0, y=t,kind='kde',ax=axes[i % 2])
#	sns.displot(datap,x=pvar0, y=t,kind='kde',cut=0,cmap='mako',fill=True,ax=axs[i % 2]) #

#	norm = plt.Normalize(datap['errsum'].max(),datap['errsum'].min())
#	sm = plt.cm.ScalarMappable(cmap='mako_r', norm=norm) #
#	sm.set_array([])
#	plt.colorbar(sm)

#	[X, Y] = np.meshgrid(data[[pvar0]],data[[t]])
#	xs=np.array(data[[pvar0]])
#	cols = np.array(len(xs))
#	X=np.array(data[[pvar0]]).reshape(0,cols)
#	Y=np.array(data[[t]]).reshape(0,cols)
#	Z=np.array(data[['errsum']]) #.values.reshape(len(X),len(Y)).T
#	Z=Z*Z
#	print (Z.shape)


#	cols = 50 #np.array(len(xs)) #.shape[0] #np.unique(xs)
#	X = xs.reshape(0, cols)
#	Y = np.array(data[[t]]).reshape(0, cols)
#	Z = np.array(data[['errsum']]).reshape(0, cols)
#	ax.tricontour(X,Y,Z)
	

#for i in range (0,3):
#	for j in range(0,2):
#		print ([i,j])#(i,j)
		#axs[i,j].sns.displot(datatoplot,x=pvar0, y=pvar1,kind='kde',cut=0,cmap='mako',fill=True)



##input for plotting, simple plotting and check color scale with scatter #####################3
plt.close('all')
#print (data)
pvar1='v1'
datatoplot=data[[pvar0, pvar1,'errsum']] #  'msesum'
#print (datatoplot)
contour=sns.displot(datatoplot,x=pvar0, y=pvar1,kind='kde',cmap='mako_r',fill=True)#cut=0,  "v3"  ,hue="errsum",binwidth=(100,100)) 
#ax.set_xlim(0,)
#ax.set_ylim(0,)
norm = plt.Normalize(datatoplot['errsum'].max(),datatoplot['errsum'].min())
sm = plt.cm.ScalarMappable(cmap='mako_r', norm=norm) #
sm.set_array([])
plt.colorbar(sm)
sms=plt.scatter(data[pvar0],data[pvar1],c=data['errsum'], s=30, cmap='Reds') #'BuGn_r' best marked with dark green
plt.colorbar(sms) # , norm=norm
plt.show()
exit()

##try seaborn contour plot###
#print (data['v5'])
plt.close('all')
datatoplot=data[['v0','v1','errsum']] #  'msesum'
contour=sns.displot(datatoplot,x="v0", y="v1",kind='kde',cut=0,cmap='mako',fill=True)#"v3"  ,hue="errsum",binwidth=(100,100)) 
#ax.set_xlim(0,)
#ax.set_ylim(0,)
norm = plt.Normalize(datatoplot['errsum'].max(),datatoplot['errsum'].min())
sm = plt.cm.ScalarMappable(cmap='mako_r', norm=norm) #
sm.set_array([])
plt.colorbar(sm)
sms=plt.scatter(data['v0'],data['v1'],c=data['errsum'], s=30, cmap='Reds') #'BuGn_r' best marked with dark green
plt.colorbar(sms) # , norm=norm


######multiple plot with seaborn############
iris=data[['v0','v1','v2','errsum_r']]  #'v0','v3','v4','errsum'
#print (iris)
g = sns.PairGrid(iris)
g.map_upper(sns.scatterplot)
g.map_lower(sns.kdeplot)
g.map_diag(sns.kdeplot, lw=3, legend=False)

#print(data['v3'])
#plt.close('all')
#plt.hexbin(data[['v0']], data[['v3']], gridsize=15) #cmap=cmap, **kwargs


#h=sns.kdeplot(data=iris, x="v0", y="v3",fill=True)#

plt.show()
exit()


################################older code#########################ignore#############
#datadir = 'Plot/data'
data = asarray(numpy.loadtxt('param.log',delimiter=' ')) #'%s/SSC_20cmab_obs.dat'%datadir
data1=asarray(numpy.loadtxt('relativeerror.txt',delimiter=' '))
#sizedata = asarray(numpy.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))




#with open ('myoutput.txt', 'rb') as fp:
#	itemlist = pickle.load(fp)
#print(itemlist[0])
#print (data)

#X1=data[:,0]
#X1.append(data1[:,0])
#err=dict()
#for i in range(0,6): 
#	err[i]=data1[:,i]
#print (err[0])
#newstuff=sorted(err[0], key=lambda x: x[1]) #smaller to largest
print("kb",data[0:10,3],data[0:10,1])
newstuff=np.sort([data[0:10,3]])
kbv= [x for _,x in sorted(zip(data[:,3],data[:,1]))]
print (newstuff,kbv)



X1=data1[:,0] #np.concatenate((data[:,0], data1[:,0])) #breakup
Y1=data1[:,1] #np.concatenate((data[:,1], data1[:,1])) #coagulation
Z1=data1[:,2] #np.concatenate((data[:,2], data1[:,2])) #density
#E1=data1[:,8] #np.concatenate((data[:,8], data1[:,8])) #totalerror

espmc=data1[:,3] #np.concatenate((data[:,3], data1[:,3])) #error SPMC
esize=data1[:,4] #np.concatenate((data[:,4], data1[:,4]))
ess=data1[:,5] #np.concatenate((data[:,5], data1[:,5]))
absvalley=data1[:,6] #np.concatenate((data[:,6], data1[:,6]))
#valley=data1[:,7] #np.concatenate((data[:,7], data1[:,7]))


fig, ax = plt.subplots()
#ax = plt.axes(projection='3d')
p=ax.scatter(X1, Y1, Z1, c=E1, s=30, cmap='Blues_r', linewidth=1)
ax.set_xlabel('breakup')
ax.set_ylabel('coagulation')
#ax.set_zlabel('')
#fig.colorbar(p)
#plt.show()
#print (len(X1))

fig, axs = plt.subplots(3, sharex=False)
axs[0].plot(X1,espmc,'--',label='kb-espmc')
axs[0].plot(X1,esize,'--',label='kb-esize')
axs[0].plot(X1,ess,label='kb-ess')

axs[1].plot(Y1,espmc,'--',label='ka-espmc')
axs[1].plot(Y1,esize,'--',label='ka-esize')
axs[1].plot(X1,ess,label='ka-ess')

axs[2].plot(Z1,espmc,'--',label='rhop-espmc')
axs[2].plot(Z1,esize,'--',label='rhop-esize')
axs[2].plot(Z1,ess,label='rhop-ess')
plt.legend()
plt.show()

