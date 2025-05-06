"""
@author: takada
Ligand binding to protein

"""

import numpy as np
import matplotlib.pyplot as plt

npoint=100
acon=10**-5
kd=10**-6

loglig=np.zeros(npoint)
bb=np.zeros(npoint)
bbapprox=np.zeros(npoint)


for i in range(npoint):
    loglig[i]=-9.+i*6./npoint
    lig=10**loglig[i]
    a=acon/kd
    l=lig/kd
    al=1+a+l
    x=0.5*(al-np.sqrt(al**2-4*a*l))
    
    bb[i]=x*kd/acon
    bbapprox[i]=lig/(lig+kd)

plt.plot (loglig, bb)
plt.plot (loglig, bbapprox, linestyle="dashed")
plt.xlabel("log_10 ([L]_total)")
plt.ylabel("B")
plt.show()




