# -*- coding: utf-8 -*-
"""
@author: takada
"""

import numpy as np
import matplotlib.pyplot as plt
"""

Numerical analysis of Morris Lecar model.
The fourth-order Runge-Kutta method for the 1st order ordinary differencial equation
with 
vector fields and nullcleines 
"""

ntime=2000
xx=np.zeros((ntime+1,2))
x=np.zeros(2)
tt=np.zeros(ntime+1)
ff=np.zeros(2)
h=0.1

vg_ca= -1.2
s_ca= 18.0
g_ca= 4.4
vn_ca=120.

vg_k=2.
s_k=30.
g_k=8.0
vn_k=-84.

g_leak=2.0

capa=20.

v_rest=-60.
tau0=25.

xmin=-80.
xmax=80.
ymin=0.
ymax=1.


# pick up i_app for your need.
#i_app=75
i_app=150.  #oscillatory
#i_app=300


#---------------------------------------------------------



def func0(x,y):
    popen_ca=0.5*(1+np.tanh((x-vg_ca)/s_ca))
    f=1/capa*(-g_ca*popen_ca*(x-vn_ca)    \
              -g_k*y*(x-vn_k)     \
              -g_leak*(x-v_rest)+i_app)
        
    return f  

def func1(x,y):
    popen_k =0.5*(1+np.tanh((x-vg_k)/s_k))
    tau=tau0/np.cosh((x-vg_k)/(2*s_k))
    f=(-1/tau)*(y-popen_k)
    return f  

def funcca(x,y):
    popen_ca=0.5*(1+np.tanh((x-vg_ca)/s_ca))
    return popen_ca  

def func(x):
    ff[0]=func0(x[0],x[1])
    ff[1]=func1(x[0],x[1])
    return ff

#----------------------------------------------------------

xx[0,0]=-60.
xx[0,1]=0.
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
plt.ylim(xmin,xmax)
plt.plot (tt, xx[:,0],color="black")
plt.xlabel("t")
plt.ylabel("x")
plt.show()

plt.figure(figsize=(16,12))
#plt.xlim(0,3)
plt.ylim(ymin,ymax)
plt.plot (tt, xx[:,1], lw=5, color="black")
plt.xlabel("t")
plt.ylabel("y")
plt.show()

plt.figure(figsize=(16,12))
#plt.xlim(0,3)
plt.ylim(0,1)
plt.plot (tt, funcca(xx[:,0],xx[:,1]), lw=5, color="black")
plt.xlabel("t")
plt.ylabel("p_open_ca")
plt.show()

# to draw vector fields

ngrid=15
ngridl=100

xg,yg = np.meshgrid(np.linspace(xmin,xmax,ngrid),np.linspace(ymin,ymax,ngrid))

su=xmax-xmin
su2=su**2
sv=ymax-ymin
sv2=sv**2
ug = func0(xg,yg)/su/np.sqrt(func0(xg,yg)**2/su2+func1(xg,yg)**2/sv2)
vg = func1(xg,yg)/sv/np.sqrt(func0(xg,yg)**2/su2+func1(xg,yg)**2)/sv2


# for drawing nullclines

xarray=np.arange(xmin,xmax,xmax/ngridl)

popen_ca=0.5*(1+np.tanh((xarray-vg_ca)/s_ca))
temp=(-g_ca*popen_ca*(xarray-vn_ca)    \
              -g_leak*(xarray-v_rest)+i_app)
    
y_xnull=1.0/(g_k*(xarray-vn_k))*temp
y_ynull=0.5*(1+np.tanh((xarray-vg_k)/s_k))

plt.figure(figsize=(16,12))
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.quiver(xg,yg,ug,vg)
plt.plot(xarray,y_xnull,"-",lw=5,color="black")
plt.plot(xarray,y_ynull,"-",lw=5,color="black")
plt.plot(xx[:,0], xx[:,1],"--" , lw=2, color="black")
plt.show()


