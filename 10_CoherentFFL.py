# -*- coding: utf-8 -*-
"""
@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

Solve differential equation for coherent FFL by
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation

"""

aa= 1.0
bb=10
cc=1.0
xon=5
dd=0.2

ntime=500
h=0.01

xx=np.zeros((ntime+1,2))
tt=np.zeros(ntime+1)
ff=np.zeros(2)

def func(x,t):   # note that y==x[0], z==x[1]
    
    xoft=xon*(1-np.exp(-bb*t))
    ff[0]=xoft**2/(1+xoft**2)-x[0]
    ff[1]=xoft**2*x[0]**2/(cc+xoft**2)/(dd+x[0]**2)-aa*x[1]
    
    return ff

xx[0,:]=0.
tt[0]=0.0

for t in range(ntime):
    tt[t+1]=(t+1)*h
    t0=t*h
    k1 = h*func(xx[t,:], t0)
    k2 = h*func(xx[t,:]+0.5*k1, t0+0.5*h)
    k3 = h*func(xx[t,:]+0.5*k2, t0+0.5*h)
    k4 = h*func(xx[t,:]+k3, t0+h)
    xx[t+1,:]=xx[t,:]+(k1+2*k2+2*k3+k4)/6

plt.plot (tt, xx[:,0])   #  y 
plt.plot (tt, xx[:,1], linestyle="dashed")   # z 
plt.xlabel("t")
plt.ylabel("y, z(t)")
plt.show()
