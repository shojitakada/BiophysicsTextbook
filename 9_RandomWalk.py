# -*- coding: utf-8 -*-
"""

@author: takada
random walk to obtain MSD curve

"""

import numpy as np
import matplotlib.pyplot as plt


ntime=100000
deltat=100
kk=0.25
#kk=0.  # for non-diffusion motion
aa=2.   # this sets the diffusion coefficient = 1.0

""" 
#copy one of the 4 cases to outside of the comment-block 
#for normal difusion
xmax=100000   #for unconstraied case
veloc=0 # without ballistic motion
pchange=0 # for a constant-biased motion

#for subdiffusion by confined motion
xmax=15 # for confined case
veloc=0 # without ballistic motion
pchange=0 # for a constant-biased motion

#for superdiffusion via a constant-biased motion
xmax=100000   #for unconstraied case
veloc=0.2  # for a short-lived ballistic or constant-biased motion
pchange=0 # for a constant-biased motion

#for superdiffusion via a constant-biased motion
xmax=100000   #for unconstraied case
veloc=0.2  # for a short-lived ballistic or constant-biased motion
pchange=0.01  # for a short-lived ballistic motion

"""

#for normal difusion
xmax=100000   #for unconstraied case
veloc=0 # without ballistic motion
pchange=0 # for a constant-biased motion


xmin=-xmax

print("k,a, xmax, veloc, pchange=",kk,aa,xmax,veloc,pchange)

xx=np.zeros(ntime+1)
tt=np.zeros(ntime+1)
msd=np.zeros(deltat)

xx[0]=0.
tt[0]=0.
rng=np.random.default_rng()


for t in range(ntime):
    tt[t+1]=(t+1)
    r=rng.random()
    if(r<kk):
        xx[t+1]=xx[t]+aa+veloc
        if(xx[t+1]>xmax):
            xx[t+1]=xmax
    elif(r<kk*2):
        xx[t+1]=xx[t]-aa+veloc
        if(xx[t+1]<xmin):
            xx[t+1]=xmin
    else:
        xx[t+1]=xx[t]+veloc
    
    r2=rng.random()
    if(r2<pchange):     #change the velocity with this probability
        veloc=-veloc

#tt=np.log(tt)
#plt.xlim([-0.0001,0.005])
plt.plot (tt, xx)
plt.xlim(0,1000)
plt.ylim(-40,40)
plt.xlabel("t")
plt.ylabel("x")
plt.show()

for t in range(deltat):
    sum=0
    for t2 in range(ntime-deltat):
        sum+=(xx[t2]-xx[t2+t])**2
    msd[t]=sum/(ntime-deltat)


plt.plot (tt[:deltat], msd)
plt.ylim(-10,200)
plt.xlabel("time")
plt.ylabel("MSD")
plt.show()

        

