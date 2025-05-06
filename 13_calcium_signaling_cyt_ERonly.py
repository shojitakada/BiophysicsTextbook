# -*- coding: utf-8 -*-
"""
@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

Solve the calcium signaling via ER membrane
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation

"""

ntime=3000
xx=np.zeros((ntime+1,2))
x=np.zeros(2)
tt=np.zeros(ntime+1)
ff=np.zeros(2)
h=0.01

fi=0.01   #based on Fall
Vi=4.0   #based on Fall
Leak=0.37   #based on Fall
Ki=1.0   #based on Fall
Ka=0.4 /2   #modified 

#P_IP3R=100000*0.95   #modified 
P_IP3R=100000*0.9   #modified 

#caT=1.7    # modified
caT=2.0     #based on Fall
    
Vserca=400.   #based on Fall
Kserca=0.2   #based on Fall

AA=0.5 *2.9    #modified
Kd=0.4 /2   #modified 
    
IP3=0.8
sigma=0.185   #based on Fall


def func(x):


    cai=x[0]
    h=x[1]
    caER= (caT-cai)/sigma

    fact_IP3=IP3/(IP3+Ki)
    fact_ca_act=cai/(cai+Ka)
    Jout_ER=(Leak + P_IP3R * (fact_IP3 * fact_ca_act *h)**4 )*(caER-cai)     
    Jin_ER=Vserca*cai**2/(cai**2+Kserca**2)
    ff[0]=fi/Vi *(Jout_ER-Jin_ER)
#    ff[1]=AA*(Kd-(cai+Kd)*h)
    ff[1]=AA*(Kd/(cai+Kd)-h)
    
    return ff

#------------------------------------------------------------------------------
xx[0,0]=0.047
xx[0,1]=0.816
tt[0]=0.0

for t in range(ntime):
    tt[t+1]=t*h
    k1 = h*func(xx[t,:])
    k2 = h*func(xx[t,:]+0.5*k1)
    k3 = h*func(xx[t,:]+0.5*k2)
    k4 = h*func(xx[t,:]+k3)
    xx[t+1,:]=xx[t,:]+(k1+2*k2+2*k3+k4)/6


print("xx=",xx[ntime-1,:])

plt.plot (tt, xx[:,0])
plt.xlabel("time (s)")
plt.ylabel("[Ca2+]cyt(t)")
plt.show()

plt.plot (tt, xx[:,1])
plt.xlabel("time (s)")
plt.ylabel("P_IP3R_d_not(t)")
plt.show()
