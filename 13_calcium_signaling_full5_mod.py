# -*- coding: utf-8 -*-
"""
@author: takada
"""

import numpy as np
import math
import matplotlib.pyplot as plt
"""

Solve the full calcium signaling with plasma and ER membrane. 
5 variables are coupled. 
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation

"""

ntime=3000*1
xx=np.zeros((ntime+1,5))
x=np.zeros(5)
tt=np.zeros(ntime+1)
ff=np.zeros(5)
h=0.01



fi=0.01   #based on Fall
Vi=4.0   #based on Fall
Leak=0.37   #based on Fall
Ki=1.0   #based on Fall
Ka=0.4 /2   #modified

P_IP3R=100000*0.95   #modified 
P_IP3R=100000*0.9   #modified 
    
Vserca=400.   #based on Fall 
Kserca=0.2    #based on Fall

AA=0.5*2.9  #modified
Kd=0.4 /2   #modified  
    
IP3=0.8 
#IP3=0.6   # a trial

V_PMCA=400.    #based on Fall
K_PMCA=0.3

# parameters for modified Morris-Lecar    
VK  =-85.   #based on Fall
Vca =120.   #based on Fall
KKca=0.5   #based on Fall

gK  =20.   #based on Fall
gca =20.    #based on Fall
gKca=8.   #based on Fall
#gKca=0.   # this is to decouple the Morris-Lecar part from others

phi =12.0   #based on Fall

v1  =-3.   #based on Fall
v2  = 30.   #based on Fall
v3  =-20.   #based on Fall
v4  = 30.   #based on Fall

# conversion parameter
alpha =0.2     #based on Fall

# geometric ratios
sigma= 0.185    # ratio of volumes based on Fall
epsilon = 0.01   # ratio of surfaces (based on Fall)
epsilon = 0.01/100   # modified

def func(x):
   
    volt=x[0]
    ww  =x[1]
    cai =x[2]
    hh  =x[3]
    caT =x[4]
    
    caER= (caT-cai)/sigma
    fact_IP3=IP3/(IP3+Ki)
    fact_ca_act=cai/(cai+Ka)
    Jout_ER=(Leak + P_IP3R * (fact_IP3 * fact_ca_act *hh)**4 )*(caER-cai)    
    Jin_ER=Vserca*cai**2/(cai**2+Kserca**2)
    
    
    minfi=0.5*(1+math.tanh((volt-v1)/v2))
    winfi=0.5*(1+math.tanh((volt-v3)/v4))
    tau=1./math.cosh((volt-v3)/(2*v4))

    Jout= V_PMCA * cai**2/(K_PMCA**2 + cai**2)
    Jin= -alpha*gca*minfi*(volt-Vca)
    
    ff[0]=-gca*minfi*(volt-Vca)-(gK*ww+ gKca*cai**4/(cai**4+KKca**4))*(volt-VK)
    ff[1]=phi*(winfi-ww)/tau    
    ff[2]=fi/Vi *(Jout_ER-Jin_ER) + epsilon* fi/Vi* (Jin-Jout)
    ff[3]=AA*(Kd/(cai+Kd)-hh)
    ff[4]=epsilon * fi/Vi *(Jin-Jout)
    
    return ff

#------------------------------------------------------------------------------
xx[0,0]=-31.1   # Volt
xx[0,1]=0.5414     # ww
xx[0,2]=0.0533     # cai
xx[0,3]=0.6683     # hh
xx[0,4]=2.085    #caT

tt[0]=0.0

for t in range(ntime):
    tt[t+1]=t*h
    k1 = h*func(xx[t,:])
    k2 = h*func(xx[t,:]+0.5*k1)
    k3 = h*func(xx[t,:]+0.5*k2)
    k4 = h*func(xx[t,:]+k3)
    xx[t+1,:]=xx[t,:]+(k1+2*k2+2*k3+k4)/6
    
print("final=",xx[ntime-1,:])

plt.plot (tt, xx[:,0])
plt.xlabel("time (s)")
plt.ylabel("Volt(mV)")
plt.show()

plt.plot (tt, xx[:,2])
plt.xlabel("time (s)")
plt.ylabel("[Ca2+]cyt (micro M)")
plt.show()

plt.plot (tt, xx[:,3])
plt.xlabel("time (s)")
plt.ylabel("P_IP3R_d_not(t)")
plt.show()

plt.plot (tt, (xx[:,4]-xx[:,2])/sigma)
plt.xlabel("time (s)")
plt.ylabel("[Ca2+]ER (micro M)")
plt.show()

plt.plot (tt, xx[:,4])
#plt.ylim(2.0,2.01)
plt.xlabel("time (s)")
plt.ylabel("[Ca2+]total (micro M)")
plt.show()