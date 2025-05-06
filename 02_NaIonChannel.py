# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 17:49:01 2021

@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

Na Ion Channel in three state model
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation

"""

ntime=1000
xx=np.zeros((ntime+1,2))
x=np.zeros(2)
tt=np.zeros(ntime+1)
ff=np.zeros(2)
h=0.01/2


def func(x):

#Folding: fast initial step, followed by slow final step
    k1=4.0
    k1min=0.11
    k2=0.97
    k2min=0.0027

    ff[0]=-k1*x[0]+k1min*x[1]
    ff[1]=k1*x[0]-k1min*x[1]-k2*x[1]+k2min*(1-x[0]-x[1])

    return ff

xx[0,0]=1.
xx[0,1]=0.
tt[0]=0.0

for t in range(ntime):
    tt[t+1]=(t+1)*h
    k1 = h*func(xx[t,:])
    k2 = h*func(xx[t,:]+0.5*k1)
    k3 = h*func(xx[t,:]+0.5*k2)
    k4 = h*func(xx[t,:]+k3)
    xx[t+1,:]=xx[t,:]+(k1+2*k2+2*k3+k4)/6


#plt.xlim([-0.1,5.])
plt.plot (tt, xx[:,0])
plt.plot (tt, xx[:,1], linestyle="dashed")
plt.plot (tt, 1-xx[:,0]-xx[:,1], linestyle="dotted")
plt.xlabel("t")
plt.ylabel("probability")
plt.show()

#"""
#tt=np.log(tt)
#plt.xlim([-1,20.])
plt.ylim([-2,0.])
plt.plot (tt, np.log(xx[:,0]))
plt.plot (tt, np.log(xx[:,1]), linestyle="dashed")
plt.plot (tt, np.log(1-xx[:,0]-xx[:,1]), linestyle="dotted")
plt.xlabel("t")
plt.ylabel("probability")
plt.show()
#"""
