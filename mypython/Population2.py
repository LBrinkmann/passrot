#! /usr/bin/env python

import sys
import os
import argparse
import scipy.ndimage
import scipy.io
import scipy.sparse
import numpy as np 
import copy
import shutil
from subprocess import call, Popen
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
from scipy.optimize import fmin,minimize
from math import exp,sqrt
import random




# initial conditions
y0 = [0.0,0.11,0.0,0.0]
#start_time,end_time,intervals = -1000,1000000,0.1

start_time,end_time,intervals = 0,1000000,1


exptimes=np.array([100,178,316,562,1000,1780,3160,5620,10000,17800,31600,56200,100000,178000,316000,1000000])
#exptimes=np.array([100,178,316,562,1000])

plottimes=np.array(range(0,1000000,1))

dT = 147
#dT = 80
T0 = 1
Nph = 0.11

sig=dT/2.355

k_r0_r1 = 0.00237
k_r1_r0 = 0.0007
k_r1_r2 = 0.0000606
        

def eq(par,initial_cond,start_t,end_t,incr):
    #-time-grid-----------------------------------
    t  = np.arange(start_t, end_t+incr,incr)
    #differential-eq-system----------------------
    def funct(y,t):
        k_r0_r1,k_r1_r0,k_r1_r2,sig,T0,Nph = par
        K = np.array(\
                     [[0,0,0,0],\
                      [0,-k_r0_r1,k_r1_r0,0],\
                      [0,k_r0_r1,-k_r1_r0-k_r1_r2,0],\
                      [0,0,k_r1_r2,0]])
        return np.dot(K,y)+[0,0,0,0]
    #integrate------------------------------------
    ds = integrate.odeint(funct,initial_cond,t,h0=0.05)
    return ds,t

par = k_r0_r1,k_r1_r0,k_r1_r2,sig,T0,Nph

ds,t = eq(par,y0,start_time,end_time,intervals)

t2 = np.insert(t,0,np.arange(-200,0,1))
tzeros = np.zeros(4*200)
tzeros = tzeros.reshape((-1,4))
ds2 = np.insert(ds,0,tzeros,axis=0)
ds = scipy.ndimage.filters.gaussian_filter(ds2,[147/2.3,0])
#ds = ds2
t = t2


timeindex = np.searchsorted(t, exptimes)
texp = t[timeindex]
dsexp = ds[timeindex,:]

print t
print ds

fig = plt.figure()
ax = fig.add_subplot(2,1,1)

plt.plot(t[::1],ds[::1,0],'k',t[::1],ds[::1,1],'r',t[::1],ds[::1,2],'b',t[::1],ds[::1,3],'g',texp,dsexp[:,1],'.r',texp,dsexp[:,2],'.b',texp,dsexp[:,3],'.g')
ax.set_xscale('log')
plt.title('Occupations [log]')
ax.set_xlim([1,1000000])
ax = fig.add_subplot(2,1,2)
plt.plot(t[::1],ds[::1,0],'k',t[::1],ds[::1,1],'r',t[::1],ds[::1,2],'b',t[::1],ds[::1,3],'g',texp,dsexp[:,1],'.r',texp,dsexp[:,2],'.b',texp,dsexp[:,3],'.g')
ax.set_xscale('linear')
ax.set_xlim([-100,1000])
ax.set_ylim([0,0.15])
plt.title('Occupations [lin]')
plt.savefig('Occupations.png')
 
