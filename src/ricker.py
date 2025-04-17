import matplotlib.pyplot as plt
import numpy as np
import sys,os

def wavelet(fm,t,t0):
    t=np.pi*fm*(t-t0)
    t=t**2
    return (1-2*t)*np.exp(-t)


fm=20
t0=0.1
t=np.arange(0,1,0.001)
wt=wavelet(fm,t,t0)
plt.plot(t,wt,linewidth=1,color	='black')
plt.fill_between(t,0,wt,where=wt>=0,facecolor='black')
plt.xlim((0,1))
plt.savefig('ricker.svg',dpi=300,format='svg')
plt.show()

