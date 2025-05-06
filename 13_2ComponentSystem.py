# -*- coding: utf-8 -*-
"""
@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

Solve two component system by
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation

"""

Lig    = 1.
pulseduration = 2.

Atotal = 6.0
Rtotal = 1.0
kk1     = 0.1
kk1m    = 0.01
kk2     = 1.5
kk3     = 0.05



ntime=10000
h=0.01

xx=np.zeros((ntime+1,2))
tt=np.zeros(ntime+1)
ff=np.zeros(2)
x=np.zeros(2)

ligand=np.zeros(ntime+1)

def func(x,time):   # note that y==x[0], z==x[1]
    
    if time > 0. and time < pulseduration:
        L=Lig
    else:
        L=0.
    
    ff[0]=kk1*L*(Rtotal-x[0])-kk1m*x[0]
    ff[1]=kk2*(Atotal - x[1])*x[0] - kk3*x[1]
    
    return ff

xx[0,:]=0.
tt[0]=0.0

for t in range(ntime):
    tt[t+1]=(t+1)*h
    t0=t*h

    if tt[t+1]> 0. and tt[t+1]< pulseduration:
        ligand[t+1]=Lig


    k1 = h*func(xx[t,:], t0)
    k2 = h*func(xx[t,:]+0.5*k1, t0+0.5*h)
    k3 = h*func(xx[t,:]+0.5*k2, t0+0.5*h)
    k4 = h*func(xx[t,:]+k3, t0+h)
    xx[t+1,:]=xx[t,:]+(k1+2*k2+2*k3+k4)/6


plt.plot (tt, xx[:,0]/Rtotal)   #  concentration of RL
plt.plot (tt, xx[:,1]/Atotal, linestyle="dashed")   # concentration of A* 
plt.plot (tt, ligand[:]/Lig, linestyle="dotted")
plt.xlabel("t")
plt.ylabel("y, z(t)")
plt.show()
