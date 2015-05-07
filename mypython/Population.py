#! /usr/bin/env python

import sys
import os
import argparse
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

        

def eq_anfin(par,initial_cond,start_t,end_t,incr):
    #-time-grid-----------------------------------
    t  = np.arange(start_t, end_t+incr,incr)
    #differential-eq-system----------------------
    def funct(y,t):
        k_r0_r1,k_r1_r0,k_r1_r2,k_r2_b,sig,T0,Nph = par
        ph = -(t-T0)/(sig**2)*Nph*exp(-0.5*((t-T0)**2)/(sig**2))
        ph2 = Nph*exp(-0.5*((t-T0)**2)/(sig**2))
        a = 1.0/(sig*sqrt(2*np.pi))
        K = np.array(\
                     [[0,0,0,0,0],\
                      [0,-k_r0_r1,k_r1_r0,0,0],\
                      [0,k_r0_r1,-k_r1_r0-k_r1_r2,0,0],\
                      [0,0,k_r1_r2,-k_r2_b,0],\
                      [0,0,0,k_r2_b,0]])
        return np.dot(K,y)+[ph,a*ph2,0,0,0]
    #integrate------------------------------------
    ds = integrate.odeint(funct,initial_cond,t,h0=0.05)
    return ds,t


def eq_ihee(par,initial_cond,start_t,end_t,incr):
    #-time-grid-----------------------------------
    t  = np.arange(start_t, end_t+incr,incr)
    #differential-eq-system----------------------
    def funct(y,t):
        k_it_pr1,k_it_ict,k_ict_pr2,k_pr1_pb,k_pr2_pb,sig,T0,Nph = par
        ph = -(t-T0)/(sig**2)*Nph*exp(-0.5*((t-T0)**2)/(sig**2))
        ph2 = Nph*exp(-0.5*((t-T0)**2)/(sig**2))
        a = 1.0/(sig*sqrt(2*np.pi))
        K = np.array(\
                     [[0,0,0,0,0,0],\
                      [0,-k_it_pr1-k_it_ict,0,0,0,0],\
                      [0,k_it_ict,-k_ict_pr2,0,0,0],\
                      [0,k_it_pr1,0,-k_pr1_pb,0,0],\
                      [0,0,k_ict_pr2,0,-k_pr2_pb,0],\
                      [0,0,0,k_pr1_pb,k_pr2_pb,0]])                 
        return np.dot(K,y)+[ph,a*ph2,0,0,0,0]
    #integrate------------------------------------
    ds = integrate.odeint(funct,initial_cond,t,h0=0.05)
    return ds,t




ratefuncdict = {
  'anfinrud': (eq_anfin,4),
  'ihee': (eq_ihee,5)
}
    


def get_populations(output,par,times,model):
    
    end_time = times[-1]
    start_time = -1000
    intervals = 0.1
    thisratefunction = ratefuncdict[model]
    ratefunction = thisratefunction[0]
    npar = thisratefunction[1]
    y0 = [0.0]*npar
    ds,t = ratefunction(par,y0,start_time,end_time,intervals)
    timeindex = np.searchsorted(t,times)
    texp = t[timeindex]
    dsexp = ds[timeindex,:]

    if npar == 4:
    
        fig = plt.figure()
        ax = fig.add_subplot(2,1,1)
        plt.plot(t[::10],ds[::10,0],'k',t[::10],ds[::10,1],'r',t[::10],ds[::10,2],'b',t[::10],ds[::10,3],'g',texp,dsexp[:,1],'.r',texp,dsexp[:,2],'.b',texp,dsexp[:,3],'.g')
        ax.set_xscale('log')
        plt.title('Occupations [log]')
        ax.set_xlim([1,1000000])
        ax = fig.add_subplot(2,1,2)
        plt.plot(t[::10],ds[::10,0],'k',t[::10],ds[::10,1],'r',t[::10],ds[::10,2],'b',t[::10],ds[::10,3],'g',texp,dsexp[:,1],'.r',texp,dsexp[:,2],'.b',texp,dsexp[:,3],'.g')
        ax.set_xscale('linear')
        ax.set_xlim([-100,1000])
        ax.set_ylim([0,0.15])
        plt.title('Occupations [lin]')
        plt.savefig(output)

    elif npar == 5:
        fig = plt.figure()
        ax = fig.add_subplot(2,1,1)
        plt.plot(t[::10],ds[::10,0],'k',t[::10],ds[::10,1],'r',t[::10],ds[::10,2],'c',t[::10],ds[::10,3],'b',t[::10],ds[::10,4],'g',texp,dsexp[:,1],'.r',texp,dsexp[:,2],'.c',texp,dsexp[:,3],'.b',texp,dsexp[:,4],'.g')
        ax.set_xscale('log')
        plt.title('Occupations [log]')
        ax.set_xlim([1,1000000])
        ax = fig.add_subplot(2,1,2)
        plt.plot(t[::10],ds[::10,0],'k',t[::10],ds[::10,1],'r',t[::10],ds[::10,2],'c',t[::10],ds[::10,3],'b',t[::10],ds[::10,4],'g',texp,dsexp[:,1],'.r',texp,dsexp[:,2],'.c',texp,dsexp[:,3],'.b',texp,dsexp[:,4],'.g')
        ax.set_xscale('linear')
        ax.set_xlim([-100,1000])
        ax.set_ylim([0,0.15])
        plt.title('Occupations [lin]')
        plt.savefig(output)
        
    return texp,dsexp[:,1:]
 
