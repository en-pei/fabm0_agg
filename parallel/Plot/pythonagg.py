# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:53:44 2020

@author: LiE
"""
import numpy as np
import matplotlib.pyplot as plt
import cmath
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib import rc
#from matplotlib.colors import TwoSlopeNorm
import matplotlib.pylab as pylab
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as tick
from sys import argv
from matplotlib.ticker import LinearLocator, FormatStrFormatter,PercentFormatter
from scipy.stats import expon
from scipy.integrate import quad
#import seaborn as sns
from scipy.special import gamma, factorial
from scipy import optimize
import pickle

idxa=[]
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx#, array[idx]

#nsteps      =12*4   #steps in 12 hours    
#hours  =np.arange(0,nsteps,1)/(nsteps/12)
class pythonagg(object):    
    def __init__(self,model):
        self.C0      =0.3 #0.5    #kg m-3
        self.Dp      =10e-6
        self.fs      =np.pi/6
        self.mu      =8.90*1e-4  #dynamic viscosity of water 25 degreeC
        self.rho_p   =2650 #kg m-3
        self.rho_w   =1000
        self.h       =0.28 #1 #water depth
        self.N       =10
        self.z       =0.05 #np.arange(0,self.h,self.h/self.N)#0.15 #0.05 #np.arange(0,1,0.1) #0.2   #depth m #0.15 from water surface is the lasersizer; which 0.05 is the OBS
        self.nf      =2 #np.arange(2.1,3.1,0.1)
        self.nu      =self.mu/self.rho_w #8.9*1e-7 #kinematic viscosity of continuous phase
        self.epsilon =1e-8     #turbulent energy disspation rate 
        self.C      =55 #3 #0.3 #forced C
        #self.D       =np.arange(10,1000,1)*1e-6 #100e-6   #m #as input
        if model=='dynamical':
            self.D       =100*1e-6 #initializing
        elif model=='equi':
            self.DN      =np.arange(0,100,1)
            self.D       =(self.DN*10+10)*1e-6 #100e-6   #m #as input
        else:
            print('unknown model type for initialization')
        #other initialization for calculations
        self.E0      =9*24*3.7*1e-6#10 folded E0 for a larger dC_r  #kg m-2 s-1 =3.7*1e-7 g cm-2 s-1 Maerz 2009
        self.tau_c   =0.29*0.01  #0.36 N m-2 0.29 Feb, 0.36 July Maerz 2009
        self.tau0    =0
        self.tau     =0#np.zeros(self.N)
        self.G=0
        self.ws=0
        self.xtime=0.0
        self.C1      =0#np.zeros(self.N)
        self.dD_c    =0#np.zeros(self.N)
        self.dD_b    =0#np.zeros(self.N)
        self.dD_d    =0#np.zeros(self.N)
        self.dD_s    =0#np.zeros(self.N)
        self.dD_r    =0#np.zeros(self.N)
        self.dD      =0#np.zeros(self.N)
        self.De      =1e-4#0#np.zeros(self.N)
        self.Deidx   =0#np.zeros(self.N)
        self.dCr     =0
        self.dCs     =0
        self.dC      =0     
        self.W=0#initializing
        self.k_A=14
        self.k_B=14*1e3
        self.e_c    =0.5 #20  #45 #40 #20#70 #efficiency of coagulation
        self.e_b    =25 #28 #55 #190 #150 #200  #efficiency of breakup
        self.p      =1
        self.q      =1/2
        self.Dr     =240*1e-6 #108*1e-6   #34 108 Feb, 34 July 
        self.dDcarr=[]
        self.dDbarr=[]
        self.dDsarr=[]
        self.dDdarr=[]
        self.dDrarr=[]
        self.Dearr=[]
        self.Garr=[]
        self.dtime=[]
        self.Dnewarr=[]
        self.dDarr=[]
        self.Darr=[]
        self.dCsarr=[]
        self.dCrarr=[]
        self.Carr=[]
        self.dCarr=[]
        self.tauarr=[]
        self.k_Aarr=[]
        
#rho_a   =rho_w+165
#Delta_rho_a=(rho_p-rho_w)*(Dp/D)**(3-nf) #here 
    def calctau(self,hours):
        omega       =np.pi*2/4
        psi         =-0.5*np.pi
        self.tau0   =1/8*(0.02+0.79 + 0.79*np.sin(omega*hours + psi)) #if 0.02 to 1.6 
        #plt.plot(hours,tau0)
        self.tau    =self.tau0*(1-self.z/self.h)
        wsPr        =0.12
        Pr          =wsPr/(0.4*np.sqrt(self.tau/self.rho_w))  #von karman constant 0.35 to 0.42 #derive when tau_c, and Pr=1.2 #Rouse number 
        y0          =0.1 #reference depth of C0
        self.C1     =self.C0*(((self.h-self.z)/self.z)*(y0/(self.h-y0)))**Pr  #equation 2.32 from Breugem 2012
        #self.n      =1/(self.fs*self.rho_p)*self.Dp**(self.nf-3)*self.D**(-self.nf)*self.C 
        self.epsilon=4.76/(np.sqrt(self.z*self.h))*np.exp(-3*self.z/self.h)*(self.tau/self.rho_w)**1.5/self.h *0.01 #*0.01 to have 1.6 tau ~ G 1.(6) #Johnson and Cowen 2017 4.76 smooth 12 rough
        self.G      =np.sqrt(self.epsilon/self.nu) #*0.5 #changed the amplitude to more similar to fabm.agg
        self.Gsin   =0.02+0.79 + 0.79*np.sin(omega*hours + psi) #if 0.02 to 1.6 
        return self.G
            
    def calcdDc(self,model):
        #constant=(1/(self.fs*self.rho_p)*self.Dp**(self.nf-3)*self.C)
        #dD_ba=[]
        #dD_ra=[]
        self.k_A    =self.e_c*np.pi*self.Dp**(self.nf-3)/(2*np.sqrt(15)*self.fs*self.nf*self.rho_p)#*0.1   #14.6 winterwerp 1998 T73
        if model=='dynamical':
#            self.ntot   =self.n #????#constant*1/(1-self.nf)*(self.D)**(1-self.nf)-constant*1/(1-self.nf)*(self.D)**(1-self.nf)
            self.dD_c=self.k_A*self.C*self.G*self.D**(4-self.nf)
        elif model=='equi':
#            self.ntot   =constant*1/(1-self.nf)*(max(self.D))**(1-self.nf)-constant*1/(1-self.nf)*(min(self.D))**(1-self.nf)
            self.dD_c=self.k_A*self.C*np.outer(self.G,(self.D**(4-self.nf)))
        else:
            print('unrecognized model for coagulation')
        return self.dD_c

    def calcdDb(self,model):
        #breakup term
        a           =1  #exponential parameter for breakup term
        
        Fy          =10e-10     #yield strength of flocs, estimated in winterwerp 1998
        self.k_B    =a*self.e_b/self.nf*(self.mu/Fy)**self.q*self.Dp**(-self.p) * 0.0002968 #14.0 10*3 winterwerp 1998 T73
        self.k_Bp   =a*self.e_b*self.Dp**(-self.p)     #kB'
        #dD_bp   =-k_Bp*(mu/Fy)**q/nf*G**(q+1)*(D-Dp)**p*D**(2*q+1)  #checking if the same as dD_b
        if model=='dynamical':    
            self.dD_b   =-self.k_B*self.G**(self.q+1)*(self.D-self.Dp)**self.p*self.D**(2*self.q+1)
        elif model=='equi':
            self.dD_b   =-self.k_B*np.outer((self.G**(self.q+1)),((self.D-self.Dp)**self.p*self.D**(2*self.q+1)))    
        else:
            print('unrecognized model for breakup')
        return self.dD_b
        #for Gi in G:
            #dD_b    =-k_B*G**(q+1)*(D-Dp)**p*D**(2*q+1)  #1D equation
        #    dD_b    =-k_B*Gi**(q+1)*(D-Dp)**p*D**(2*q+1)
        #    dD_ba.append(dD_b)
        
        
    def calcdDd(self):  #differential settling term
        self.nf1         =3
        g           =9.81
        Re          =0 #1e-10
        alpha       =1
        beta        =1
        #wsp     =alpha/(beta*18*mu)*(rho_p-rho_w)*g*Dp**(6-2*nf)*D**(2*nf-4)/(1+0.15*Re**0.687) #ws to deltarho_p
        self.W      =alpha/(beta*18*self.mu)*(self.rho_p-self.rho_w)*g*self.Dp**(3-self.nf1)/(1+0.15*Re**0.687)*1
#       Dmin         =self.Dp
#       self.Di      =1/2*(self.D-Dmin)+Dmin
#       dD_d0        =-self.e_c*np.pi*self.Dp**(self.nf-3)/(2*self.fs*self.nf*self.rho_p)*self.W*(self.D**(self.nf-1)-self.D**(3-self.nf)*self.Di**(2*self.nf-4))*self.C #1D equation
#       #W1      =alpha/(beta*18*mu)*(rho_p-rho_w)*g*Dp**(6-2*nf1)/(1+0.15*Re**0.687) #sinking but nf1=3
#       E_d          =1/32*self.D**4-3/16*Dmin*self.D**3-1/16*Dmin**2*self.D**2  #expected value, poly equation of D
#       dD_d0        =(-self.e_c*np.pi*self.W/(2*self.ntot))*(1/32*self.D**4-3/16*Dmin*self.D**3-1/16*Dmin**2*self.D**2)  #based on n
#       self.dD_d         =(self.e_c*np.pi*self.W/(128*self.ntot))*(3/4*self.D**4-Dmin*self.D**3-Dmin**2*self.D**2)  #based on n
        self.dD_d     =2*np.pi*(np.e-1)/(np.e**3)*self.e_c*self.W/(self.nf*self.fs*self.rho_p)*self.D**(5-self.nf)*self.C
        return self.dD_d
       #dD_d1   =(-e_c*np.pi/2*fs*rho_p*(1/C)*1/(18*mu)*(rho_p-rho_w)*g*Dp**(3-nf)*D**nf)*(1/32*D**4-3/16*Dmin*D**3-1/16*Dmin**2*D**2) #based on C
       #dD_d2   =(-e_c*np.pi/2*fs*rho_p*(1/C)*1/(18*mu)*(rho_p-rho_w)*g*Dp**(3-nf)*D**nf)*E_d
       
       
    def calcdDs(self): #sinking velocity ws
        #ws_W    =W*D**(nf-1)
        #ws0     =alpha/(beta*18*self.mu)*(self.rho_p-self.rho_w)*g*self.Dp**(3-nf1)/(1+0.15*Re**0.687)*self.D**(nf1-1) #1D equation
        self.ws     =self.W*self.D**(self.nf1-1) #1D equation
        #for nfi in nf:    
        #    ws  =alpha/(beta*18*mu)*(rho_a-rho_w)*g*np.outer((Dp**(3-nfi)),(D**(nfi-1)))/(1+0.15*Re**0.687)
        #dD_s0    =-alpha/(beta*18*self.mu)*(self.rho_p-self.rho_w)*g*self.Dp**(6-2*self.nf)/(1+0.15*Re**0.687)*self.D**(2*self.nf-3)*1/z #1D equation
        #dD_s        =-1/self.ntot*(self.rho_p-self.rho_w)*g/(self.z*18*self.mu)*4*self.D**3
        dD_s        =(1-2/np.e)*np.log(1-2/np.e)*self.D
        return dD_s

       
    def calcdDr(self,model): #resuspension    
        #for taui in self.tau:
        if self.tau<self.tau_c:
            self.tau=self.tau_c
        #print("tau too small, no resuspension")
        else:# taui>=tau_c:
            self.tau=self.tau
        #dD_r0   =E0/z*(tau/tau_c-1)*(Dr-D)/(nf*C) #1D equations
        #dD_r1   =1/ntot*E0/z*(tau/tau_c-1)*Dp**(nf-3)/(fs*rho_p)*D**(-nf)*(Dr-D)*gamma(nf+1) #or Dp-D
        if model=='dynamical':
            #self.dD_r    =1/self.n*self.E0*(1/self.z*(self.tau/self.tau_c-1))*(self.Dp**(self.nf-3)/(self.fs*self.rho_p)*self.D**(-self.nf)*(self.Dr-self.D)*gamma(self.nf+1))
            self.dD_r    =1/self.C*self.E0*(1/self.z*(self.tau/self.tau_c-1))*(self.Dr-self.D)
        elif model=='equi':
            #self.dD_r    =1/self.n*self.E0*np.outer((1/self.z*(self.tau/self.tau_c-1)),(self.Dp**(self.nf-3)/(self.fs*self.rho_p)*self.D**(-self.nf)*(self.Dr-self.D)*gamma(self.nf+1))) #or Dp-D with gamma function with ntot
            self.dD_r    =1/self.C*self.E0*self.E0*np.outer((1/self.z*(self.tau/self.tau_c-1)),(self.Dr-self.D))
        else:
            print('unrecognized model for resuspension')
        return self.dD_r
        #dD_r    =1/n*E0/z*(tau/tau_c-1)*Dp**(nf-3)/(fs*rho_p)*D**(-nf)*(Dr-D)*gamma(nf+1) #or Dp-D with gamma function with ntot
        #dD_ra.append(dD_r)
        #dD_r2   =1/n*E0/z*(tau/tau_c-1)*Dp**(nf-3)/(fs*rho_p)*D**(-nf)*(Dr-D)

    def find_De(self):
        value = 0.0
        array=self.dD
        return find_nearest(array, value)
    
        
    def calcdCs(self):
        self.dC_s=-self.C*self.ws/self.z
        return self.dC_s
    
    def calcdCr(self):
        self.dC_r=self.E0/self.z*(self.tau/self.tau_c-1)
        return self.dC_r
        
    def runmodel(self,dt,tend,model):
        it=0
        #self.Dnew=1e-4
        #self.C=0.3#doing again
        while self.xtime<tend:
            it += 1    
            hours=self.xtime
            self.G=self.calctau(hours)
            self.dC_s=self.calcdCs()
            self.dC_r=self.calcdCr()
            self.dD_c=self.calcdDc(model)
            self.dD_b=self.calcdDb(model)
            self.dD_d=self.calcdDd()
            self.dD_s=self.calcdDs()
            self.dD_r=self.calcdDr(model)
            k=1e-3*3.2#for sinking 
            k1=0.01 #for resuspension
            k2=1 #200 #for differential settling
            self.dC=self.dC_s+self.dC_r #*1e-3
            
            if model=='equi':
                #self.C += self.dC*dt
                self.C      =3*self.G/1.41 #0.3*self.tau/0.15 #first order of C changing with tau #here
                #self.dD =self.dD_c[:]+self.dD_b[:]+self.dD_s+self.dD_d+self.dD_r[:]
                #D =(self.DN*10+10)*1e-6
                #fun=self.k_A*self.C*self.G*x**(4-self.nf)\
                #-self.k_B*self.G**(self.q+1)*(x-self.Dp)**self.p*x**(2*self.q+1)\
                #+2*np.pi*(np.e-1)/(np.e**3)*self.e_c*self.W/(self.nf*self.fs*self.rho_p)*x**(5-self.nf)*self.C\
                #+(1-2/np.e)*np.log(1-2/np.e)*x+1/self.C*self.E0*(1/self.z*(self.tau/self.tau_c-1))*(self.Dr-x)
                self.De= optimize.root(lambda x: self.k_A*self.C*self.G*x**(4-self.nf)\
                -self.k_B*self.G**(self.q+1)*(x-self.Dp)**self.p*x**(2*self.q+1)\
                +2*np.pi*(np.e-1)/(np.e**3)*self.e_c*self.W/(self.nf*self.fs*self.rho_p)*x**(5-self.nf)*self.C\
                +(1-2/np.e)*np.log(1-2/np.e)*x*k+1/self.C*self.E0*(1/self.z*(self.tau/self.tau_c-1))*(self.Dr-x),1).x
                                       ###sinking *1e-4 as it is too big
                #self.Deidx =self.find_De()
                #self.De=np.transpose(self.D[self.Deidx])
                #self.D =(self.DN*10+10)*1e-6 #rewrite D to initial array
            elif model=='dynamical':
                self.dD =self.dD_c+self.dD_b+self.dD_d*k2+self.dD_s*k+self.dD_r*k1#here #amplified differential settling
                #self.Dnew += self.dD*dt
                #self.D=self.Dnew    
                self.D += self.dD*dt            
                self.C += self.dC*dt#*1e-3 #+self.C1 #unit???????
                #self.C=self.C1
                #self.C=0.3
            self.xtime += dt/3600.0
            if it % 60 == 1:
                self.dtime.append(self.xtime)
                self.Garr.append(self.G)
                self.tauarr.append(self.tau)
                self.dDcarr.append(self.dD_c)
                self.dDbarr.append(self.dD_b)
                self.dDdarr.append(self.dD_d) 
                self.dDsarr.append(self.dD_s)
                self.dDrarr.append(self.dD_r)
                self.Dearr.append(self.De*1e6) #
                #self.Dnewarr.append(self.Dnew)
                self.dDarr.append(self.dD)
                self.Darr.append(self.D*1e6)#
                self.dCsarr.append(self.dC_s)
                self.dCrarr.append(self.dC_r)
                self.Carr.append(self.C)
                self.dCarr.append(self.dC)
                self.k_Aarr.append(self.k_A)


#%% 
d=pythonagg('dynamical')
d.runmodel(1,24,'dynamical') #til 12h, 60s time steps 
#f=open('myoutput.txt', 'w')
#print([d.dtime,d.Darr], file = f)

with open('myoutput.txt', 'wb') as fp:
	pickle.dump([d.dtime,d.Darr,d.Carr], fp)

#plt.plot(d.dtime,d.Garr)
#plt.plot(d.dtime,d.Carr) 
#%%
#plt.plot(d.dtime,d.Darr,linewidth=3) 
#plt.scatter(sdatatime,sdatasize,color='k')   
#%% check C change with time
         
#%%run with equilibrium size 
#e=pythonagg('equi')
#e.runmodel(1,24,'equi') #til 12h, 60s time steps

#%%
#d.dDcarr=np.array(d.dDcarr)
#d.dDbarr=np.array(d.dDbarr)
#d.dDdarr=np.array(d.dDdarr)
#d.dDsarr=np.array(d.dDsarr)
#d.dDrarr=np.array(d.dDrarr)
#d.dDarr=np.array(d.dDarr)
#d.Darr=np.array(d.Darr)
#d.Carr=np.array(d.Carr)
#d.dCarr=np.array(d.dCarr)
#d.dCsarr=np.array(d.dCsarr)
#d.dCrarr=np.array(d.dCrarr)
#%%
#plt.plot(d.dtime, d.dDcarr,label='coagulation')
#plt.plot(d.dtime, d.dDbarr,label='breakup')
#plt.plot(d.dtime, d.dDdarr,label='differential settling')
#plt.plot(d.dtime, d.dDsarr,label='sinking')
#plt.plot(d.dtime, d.dDrarr,label='resuspension')
#plt.plot(d.dtime, d.dDarr,label='total')
#plt.legend()
#%%
#plt.plot(d.dtime, d.dCsarr,label='sinking')
#plt.plot(d.dtime, d.dCrarr,label='resuspension')
#plt.plot(d.dtime, d.dCarr,label='total')
#plt.plot(d.dtime, d.Carr,label='C')

#test=np.array(f.Dnewarr)
#test1=test[:,0,:]         
#plt.plot(d.dtime,d.Carr)
#plt.plot(e.dtime,e.Carr) 

#%%#%% data SSC
#datadir = r'C:\Users\LiE\Documents\HZG\models\VMenpei'
#data = np.asarray(np.loadtxt('%s/SSC_20cmab_obs.dat'%datadir,delimiter=', '))
#sizedata = np.asarray(np.loadtxt('%s/D50_ESD_13cmab_obs.dat'%datadir,delimiter=', '))
#plt.plot(data[:,0]-5,data[:,2]/1000,'k^',ms=5.,label='SPM, data')#'TSM\nBale.ea2002')
#plt.plot(d.dtime,d.Carr)

#%% data size
#sdatatime=np.c_[sizedata[:,0]-5,sizedata[:,0]-5,sizedata[:,0]-5]
#sdataz=np.c_[-(0.28-sizedata[:,2]),-(0.28-sizedata[:,2]),-(0.28-sizedata[:,2])]  #-(0.28-0.13)
#sdatasize=np.c_[sizedata[:,1],sizedata[:,1],sizedata[:,1]]
#plt.plot(sdatatime,sdatasize)

#%%
fig, ax1 = plt.subplots()
color = 'tab:blue'
ax1.set_xlabel('time (h)')
ax1.set_ylabel('G (s$^{-1}$)', color=color)
ax1.plot(d.dtime,d.Garr,color=color,linestyle='-',linewidth=3,alpha=0.3)
#ax1.plot(d.dtime, d.Carr, color='purple',linestyle='-',linewidth=3) #f.Dearr,
#ax1.plot(e.dtime,e.Garr,color=color,linestyle=':',linewidth=3,alpha=0.5)
#ax1.plot(f.dtime,f.Garr,color=color,linestyle='-',linewidth=3,alpha=0.5)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'k' #tab:green
#ax2.set_ylabel('De step ($\mu$m)', color=color) #De idx for equi # we already handled the x-label with ax1
#ax2.plot(hours,idxa, color=color)
ax2.plot(d.dtime,d.Darr, color=color,linestyle='-',linewidth=3)
#ax2.plot(e.dtime,e.Dearr, color='k',linestyle='--',linewidth=3)
#ax2.plot(d.dtime, d.Carr, color='purple',linestyle='-',linewidth=3) #f.Dearr,
ax2.tick_params(axis='y', labelcolor=color)
#ax2.set_ylim(0,200)

#ax2.scatter(sdatatime,sdatasize,color='k') #data size
ax2.set_ylabel('D ($\mu$m)', color=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()



