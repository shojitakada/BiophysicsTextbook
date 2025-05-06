# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 17:49:01 2021

@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

Analysis of 2 cyclin model
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation
with 
vector fields and nullcleines

"""

ntime=20000
xx=np.zeros((ntime+1,2))
x=np.zeros(2)
tt=np.zeros(ntime+1)
ff=np.zeros(2)
h=0.1


epsilon2=0.1**2
aa= 0.1
bb= 0.1
cc= 100.
tau0= 5.0
xmax=4.
ymax=6.



#"""
aa= 0.1
bb= 1.5
cc= 1.
tau0= 5.0
xmax=4.
ymax=6.
#"""

#---------------------------------------------------------

def func0(x,y):
    f=(epsilon2+x**2)/(1+x**2)/(1+y)-aa*x
    return f  

def func1(x,y):
    f=(bb-y/(1+cc*x**2))/tau0
    return f  

def func(x):
    ff[0]=func0(x[0],x[1])
    ff[1]=func1(x[0],x[1])
    return ff

#----------------------------------------------------------

xx[0,0]=1.
xx[0,1]=1.
tt[0]=0.0

for t in range(ntime):
    tt[t+1]=t*h
    k1 = h*func(xx[t,:])
    k2 = h*func(xx[t,:]+0.5*k1)
    k3 = h*func(xx[t,:]+0.5*k2)
    k4 = h*func(xx[t,:]+k3)
    xx[t+1,:]=xx[t,:]+(k1+2*k2+2*k3+k4)/6

plt.figure(figsize=(16,12))
#plt.xlim(0,3)
plt.ylim(0,ymax)
plt.plot (tt, xx[:,0],color="black")
plt.plot (tt, xx[:,1],"--" , color="black")
plt.xlabel("t")
plt.ylabel("x, y")
plt.show()


ngrid=15
ngridl=100

xg,yg = np.meshgrid(np.linspace(0.01,xmax,ngrid),np.linspace(0.01,ymax,ngrid))

ug = func0(xg,yg)/np.sqrt(func0(xg,yg)**2+func1(xg,yg)**2)
vg = func1(xg,yg)/np.sqrt(func0(xg,yg)**2+func1(xg,yg)**2)

xarray=np.arange(0.01,xmax,xmax/ngridl)
y_xnull=(epsilon2+xarray**2)/(1+xarray**2)/(aa*xarray)-1
y_ynull=bb*(1+cc*xarray**2)

plt.figure(figsize=(16,12))
plt.xlim(0,xmax)
plt.ylim(0,ymax)
plt.quiver(xg,yg,ug,vg)
plt.plot(xarray,y_xnull,"-",lw=5,color="black")
plt.plot(xarray,y_ynull,"-",lw=5,color="black")
plt.plot(xx[:,0], xx[:,1],"--" , lw=2, color="black")
plt.show()


