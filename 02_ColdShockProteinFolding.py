# -*- coding: utf-8 -*-
"""
@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

Cold Shock protein folding in two-state model
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation

"""

ntime=50000
xx=np.zeros((ntime+1,3))
x=np.zeros(1)
tt=np.zeros(ntime+1)
ff=np.zeros(1)
h=0.001

dlist = [0, 3.3, 5.0] # list of denaturant concentrations

def func(x,denat):
    k1=565*np.exp(-2.66*denat)
    k1min=0.018*np.exp(0.419*denat)
    ff[0]=-k1*x+k1min*(1-x)
    return ff

xx[0,:]=1.
tt[0]=0.0


for dd in range(len(dlist)): 
    denat=dlist[dd]
    for t in range(ntime):
        tt[t+1]=(t+1)*h
        k1 = h*func(xx[t,dd],denat)
        k2 = h*func(xx[t,dd]+0.5*k1,denat)
        k3 = h*func(xx[t,dd]+0.5*k2,denat)
        k4 = h*func(xx[t,dd]+k3,denat)
        xx[t+1,dd]=xx[t,dd]+(k1+2*k2+2*k3+k4)/6


#plt.xlim([-0.1,5.])
plt.plot (tt, 1-xx[:,0])
plt.plot (tt, 1-xx[:,1], linestyle="dashed")
plt.plot (tt, 1-xx[:,2], linestyle="dotted")
plt.xlabel("t")
plt.ylabel("probability")
plt.show()

#log scale in x-axis
tt=np.log(tt)
#plt.xlim([-0.0001,0.005])
plt.plot (tt, 1-xx[:,0])
plt.plot (tt, 1-xx[:,1], linestyle="dashed")
plt.plot (tt, 1-xx[:,2], linestyle="dotted")
plt.xlabel("t")
plt.ylabel("probability")
plt.show()


