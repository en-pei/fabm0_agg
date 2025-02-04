import matplotlib
from matplotlib import rc
from matplotlib import ticker, cm
import netCDF4
from pylab import *
import sys
import cftime #netcdftime
import numpy
import pickle
from mpl_toolkits import mplot3d
#from astropy.visualization import simple_norm

rc('font',**{'family':'serif','serif':['FreeSans']})


#datadir = 'Plot/data'
data = asarray(numpy.loadtxt('27relativeerror.txt',delimiter=' ')) #'%s/SSC_20cmab_obs.dat'%datadir
data1=asarray(numpy.loadtxt('relativeerror.txt',delimiter=' '))
#sizedata = asarray(numpy.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))

#with open ('myoutput.txt', 'rb') as fp:
#	itemlist = pickle.load(fp)
#print(itemlist[0])
#print (data)

#X1=data[:,0]
#X1.append(data1[:,0])
X1=data1[:,0] #np.concatenate((data[:,0], data1[:,0])) #breakup
Y1=data1[:,1] #np.concatenate((data[:,1], data1[:,1])) #coagulation
Z1=data1[:,2] #np.concatenate((data[:,2], data1[:,2])) #density
E1=data1[:,8] #np.concatenate((data[:,8], data1[:,8])) #totalerror

espmc=data1[:,3] #np.concatenate((data[:,3], data1[:,3])) #error SPMC
esize=data1[:,4] #np.concatenate((data[:,4], data1[:,4]))
ess=data1[:,5] #np.concatenate((data[:,5], data1[:,5]))
absvalley=data1[:,6] #np.concatenate((data[:,6], data1[:,6]))
valley=data1[:,7] #np.concatenate((data[:,7], data1[:,7]))


fig, ax = plt.subplots()
ax = plt.axes(projection='3d')
p=ax.scatter(X1, Y1, Z1, c=E1, s=E1*30, cmap='Blues', linewidth=1)
ax.set_xlabel('breakup')
ax.set_ylabel('coagulation')
ax.set_zlabel('density')
fig.colorbar(p)
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

