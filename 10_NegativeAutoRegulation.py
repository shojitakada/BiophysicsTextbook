# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 17:49:01 2021

@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

solve differential equation of a negative auto regulation.
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation

"""

aa = 0.1
bb = 0.2

ntime=500
h=0.1


tt=np.zeros(ntime+1)
xx=np.zeros((ntime+1,2))
xx2=np.zeros((ntime+1,2))

def f(x):   # negtive auto regulation
    return 1/(1+x*x)-aa*x

def f2(x):  # without regulation
    return bb-aa*x

xx[0] = 0.
xx2[0] = 0.
tt[0] = 0.

for t in range(ntime):
    
    tt[t+1]=tt[t]+h
    k1 = h*f(xx[t])
    k2 = h*f(xx[t]+0.5*k1)
    k3 = h*f(xx[t]+0.5*k2)
    k4 = h*f(xx[t]+k3)
    xx[t+1] = xx[t]+(k1+2*k2+2*k3+k4)/6

    k1 = h*f2(xx2[t])
    k2 = h*f2(xx2[t]+0.5*k1)
    k3 = h*f2(xx2[t]+0.5*k2)
    k4 = h*f2(xx2[t]+k3)
    xx2[t+1] = xx2[t]+(k1+2*k2+2*k3+k4)/6

plt.plot (tt, xx)
plt.plot (tt, xx2, linestyle="dashed")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.show()