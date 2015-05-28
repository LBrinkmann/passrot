#! /usr/bin/env python

import sys
import argparse
from argparse import RawTextHelpFormatter
from math import sqrt,sin,cos,acos,atan2,pi,asin,acos,isnan,exp
#from subprocess import check_call,call,Popen
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from mypython.WAXSPrint import print_mapping
from mypython.ExperimentRef import ExperimentRef
from mypython.ProteinRef import ProteinRef,check_consistency
from mypython.Population import get_populations

# in eV*nm
hc = 1239.842


#exptimes=np.array([100,1000,10000,100000,1000000])


#parameter anfinrud (Schotte et. al 2012)
#files anfinrud: pR0,pR1,pR2,pB
#files jung nomencalutur: It,Icp,pR2,pB

dT = 147
T0 = 1
Nph = 0.11

sig=dT/2.355

k_r0_r1 = 0.00237      
k_r1_r0 = 0.0007       
k_r1_r2 = 0.0000606    
k_r2_b = 0.00000000122
k_b_pg = 3.84 * 10**-12

par_anfinrud = k_r0_r1,k_r1_r0,k_r1_r2,k_r2_b,k_b_pg,sig,T0,Nph

#parameter jung (Jung et al. 2013)
#files ihee: It,Ict,pR1,pR2
#APS
k_it_pr1 = 1.0/3000 
k_it_ict = 1.0/1700   
k_ict_pr2 = 1.0/20000

par_jung_aps = k_it_pr1,k_it_ict,k_ict_pr2,sig,T0,Nph

#ESRF
k_it_pr1 = 1.0/700
k_it_ict = 1.0/400
k_ict_pr2 = 1.0/6000

par_jung_esrf = k_it_pr1,k_it_ict,k_ict_pr2,sig,T0,Nph

#parameter ihee (Ihee et al. 2005)
#files ihee: Icp,pRe46q,pRcw,pB1,pB2
#files jung: Icp,pR1,pR2,pB1,pB2
k1 = 4.8 * 10**-5  
k2 = 3.7 * 10**-5
k3 = 3.0 * 10**-9
k4 = 3.3 * 10**-8
k5 = 55 * 10**-12
k6 = 100 * 10**-12
k7 = 7.1 * 10**-12

par_ihee = k1,k2,k3,k4,k5,k6,k7,sig,T0,Nph

par_dir = {
  'anfinrud': par_anfinrud,
  'jung_aps': par_jung_aps,
  'jung_esrf': par_jung_esrf,
  'ihee': par_ihee
}



            
def doublecos2weight(qin,qout):
    raise NameError('doublecos2weight: Still to be done.')


def cos2_matrix(allqin,allqabsin,allqout,allqabsout):
    difffactor=exp(-(args.t*args.D*6))
    output = np.zeros([len(allqabsin.abs),args.phi,allqabsin.nq[-1]])
    for k,(qabs,i_in,nq_in,i_out,nq_out) in enumerate(zip(allqabsin.abs,allqabsin.index,allqabsin.nq,allqabsout.i_out,allqabsout.nq_out)):
        print "Calculation for qabs = "+str(qabs)
#        if k == 1:
        if qabs != 0:
            for i,q_out in enumerate(allqout[i_out:i_out+nq_out]):
                lq=np.dot(q_out,[0,1,0])/qabs
                for j,q_in in enumerate(allqin[i_in:i_in+nq_in]):
                    mq=np.dot(q_in,args.m[0:3])/qabs
                    if args.mtype == 'linear':
                        output[k][i][j] = ((1.5*lq*lq-0.5)*(3*mq*mq-1)+1)
                    elif args.mtype == 'planar':
                        output[k][i][j] = -(1.5*lq*lq-0.5)*(1.5*mq*mq-0.5)+1
        else:
            output[k][0][0] = 1
        #if args.illustrate !=0 and (k-1)%args.illustrate == 0 and k > 0:
        if args.illustrate !=0 and k == 1:
            print_mapping(args.o+'_l2',allqin,allqabsin,allqout,allqabsout,output,k,args.m)
            #better_sphere(allqin,allqabsin,allqout,allqabsout)
            #write_sphere(allqin,allqabsin,allqout,allqabsout,output)
    return output

def dcos2_matrix(allqin,allqabsin,allqout,allqabsout):
    output = np.zeros([len(allqabsin.abs),args.phi,allqabsin.nq[-1]])
    m=np.array(args.m[0:3])
    n=np.array(args.m[3:6])
    print m, n
    if args.mtype == 'planar':
        raise NameError('Planar transition moment not supported for two transition moments.')
    for k,(qabs,i_in,nq_in,i_out,nq_out) in enumerate(zip(allqabsin.abs,allqabsin.index,allqabsin.nq,allqabsout.i_out,allqabsout.nq_out)):
        
        print "Calculation for qabs = "+str(qabs)
#        if k == 1:
        if qabs != 0:
            for i,q_out in enumerate(allqout[i_out:i_out+nq_out]):
                lq=np.dot(q_out,[0,1,0])/qabs
                lq2=lq*lq
                lt2=(1-lq2)
                for j,q_in in enumerate(allqin[i_in:i_in+nq_in]):
                    mq=np.dot(q_in,m)/qabs
                    mq2=mq*mq
                    mqvec=mq*q_in/qabs
                    mt2=1-mq2
                    mtvec=m-mqvec
                    nq=np.dot(q_in,n)/qabs
                    nq2=nq*nq
                    nqvec=nq*q_in/qabs
                    nt2=1-nq2
                    ntvec=n-nqvec
                    mtdotnt = np.dot(mtvec,ntvec)
                    output[k][i][j] = lq2*lq2*mq2*nq2+0.5*lq2*lt2*(mq2*nt2+nq2*mt2+4*mq*nq*mtdotnt)+0.25*lt2*lt2*mtdotnt*mtdotnt
                        
        else:
            output[k][0][0] = 1
        #if args.illustrate !=0 and (k-1)%args.illustrate == 0 and k > 0:
        if args.illustrate !=0 and k == 1:
            print_mapping(args.o+'_l2',allqin,allqabsin,allqout,allqabsout,output,k,args.m)
            #better_sphere(allqin,allqabsin,allqout,allqabsout)
            #write_sphere(allqin,allqabsin,allqout,allqabsout,output)

    return output

funcdict = {
  'cos2': cos2_matrix,
  'dcos2': dcos2_matrix
}




    
parser = argparse.ArgumentParser(description='Anisotropic pattern by passive rotation.\n\n'
                                 'The output files correspond to a Laser excitation perpendicular\n'
                                 'to the X-ray beam for linear polarisation. We assume a vertical\n'
                                 'laser beam orientation.\n'
                                 'In the case of linear polarisation the X-ray beam is perpendic-\n'
                                 'ular to the laser polarisation.\n'
                                 'In the case of circular polarisation the normal of the laser\n'
                                 'exciation is perpendicular to X-ray beam.\n'
                                 'In both cases the parallel case is writen into *_par.xvg.\n'
                                 , formatter_class=RawTextHelpFormatter)
parser.add_argument('-fa', nargs='+', metavar='_averagA.xvg',help='Average scattering amplitude of System A.')
parser.add_argument('-fb', nargs='+', metavar='_averagB.xvg',help='Average scattering amplitude of System B.')
parser.add_argument('-env', metavar='',help='Fourier transform of the envelope.')
parser.add_argument('-env_ab', action='store_true', help='Assume envelope files for input -fa and -fb.')
parser.add_argument('-o', metavar='.xvg',help='Output.')
parser.add_argument('-vacuum', action='store_true', help='Do vacuum calculation.')
parser.add_argument('-phi', type=int, default=100, help='Number of azimutal data points in 2D output.')
parser.add_argument('-weight', choices=['cos2','dcos2'], default='cos2', help='Weight of different orientations.')
parser.add_argument('-popfunc', choices=['ihee','anfinrud','jung_aps','jung_esrf'], default='anfinrud', help='Function for occupation.')
parser.add_argument('-m', type=float, default=[1,0,0], nargs='+', help='Vector of excitation moment(s).')
parser.add_argument('-mtype', choices=['linear','planar'], default='linear',
                   help='Type of the excitation moment.')
parser.add_argument('-polar', choices=['linear','circular'], default='linear',
                   help='Type of polarisation.')
parser.add_argument('-beam', type=float, default=12.0, help='X-ray beam energy. [keV]')
parser.add_argument('-bulk', type=float, default='nan', help='Set bulk water density.')
parser.add_argument('-D', type=float, default=0.000027, help='Rotational diffusion constant.')
parser.add_argument('-exptimes', type=int, nargs='+', default=[100,178,316], help='Exptimes.')
parser.add_argument('-numeric', action='store_true',help='Do numeric radial average.')
parser.add_argument('-legend', help='Name in legend for plot.')
parser.add_argument('-azimutal', action='store_true', help='Use azimutal average instead of horizontal and vertical cuts.')
#parser.add_argument('-model', type=float, help='Assume a homogenous protein in homogenous solvent. Value correspond to protein electron density.')
parser.add_argument('-scale', type=float, default=1, help="Scaling of the electron density.")
parser.add_argument('-illustrate', default=0, type=int, help="Illustrate the rotational averaging. 0:do not, >0: illustrate:nq%")
parser.add_argument('-plotdesign', choices=['cho','kim'], default='cho', help="Design of the plots.")


args = parser.parse_args()


exptimes = np.array(args.exptimes)

nparts=len(args.fb)
if len(args.fa)%len(args.fb) !=0:
    raise NameError('A multiple of files have to be given for A then for B.')
nspecs=len(args.fa)/len(args.fb)



texp,popexp = get_populations(args.o+"populations.png",par_dir[args.popfunc],exptimes,args.popfunc)


allA = [ProteinRef() for fa in args.fa]
for fname,A in zip(args.fa,allA):
    if not args.env_ab:
        A.read_file(fname)
    else:
        A.read_envelope(fname)
    check_consistency(allA[0],A)
    
allB = [ProteinRef() for fb in args.fb]
for fname,B in zip(args.fb,allB):
    if not args.env_ab:
        B.read_file(fname)
    else:
        B.read_envelope(fname)
    check_consistency(allA[0],B)

if args.env:
     E = ProteinRef()
     E.read_envelope(args.env)
     check_consistency(allA[0],E)



someargs = args.phi,args.beam,args.bulk,args.numeric,args.m,args.illustrate,args.D,args.polar,args.weight,args.mtype,args.legend,args.azimutal,args.scale,args.vacuum,nparts
allI = [ExperimentRef(allB[0].nqabs,allB[0].qabs,someargs) for i in range(nspecs)]

for i,I in enumerate(allI):
    if args.env:
        I.calc_q_in_intensity_error(allA[i*nparts:(i+1)*nparts],allB,E,args.o)
    else:
        I.calc_q_in_intensity_error(allA[i*nparts:(i+1)*nparts],allB,None,args.o)

allIout = [ExperimentRef(allB[0].nqabs,allB[0].qabs,someargs) for t in exptimes]
pp = PdfPages(args.o+"_timeline.pdf")
fig = plt.figure(figsize=(7, 6))
ax_avg = fig.add_subplot(121)
ax_avg.set_title("radial averaged",fontsize = 15)
if args.plotdesign == 'cho': 
    ax_HminV = fig.add_subplot(122)
    ax_HminV.set_title("horizontal - vertical",fontsize = 15)
else:
    ax_HminV = None

shift_avg=None
shift_HminV=None
for i,(expt,Iout) in enumerate(zip(exptimes,allIout)):
    Iout.combine_intensities(allI,popexp[i,:])
    Iout.init_q()
    if not args.numeric:
        Iout.do_passiv_rot(expt)
        Iout.write_intensity(expt,args.o+"_"+str(expt))
        shift_avg,shift_HminV = Iout.add_plots(ax_avg,ax_HminV,i,expt,shift_avg,shift_HminV,str(expt)+"ps",args.plotdesign)
    else:
        Iout.do_passiv_rot_numeric(funcdict[args.weight])
        Iout.write_numeric(args.o+"_"+str(expt))

if args.plotdesign == 'kim': 
    ax_avg.get_yaxis().set_visible(False)
    ax_avg.set_xscale('log')
    ax_avg.set_xlim(0.5,allB[0].qabs.abs[-1])

if args.plotdesign == 'cho': 
    ax_avg.get_yaxis().set_visible(False)
    #ax_avg.autoscale(axis='y')
    ax_avg.set_xlim(allB[0].qabs.abs[0],allB[0].qabs.abs[-1])
    ax_HminV.get_yaxis().set_visible(False)
    #ax_HminV.autoscale(axis='y')
    ax_HminV.set_xlim(allB[0].qabs.abs[0],allB[0].qabs.abs[-1])

pp.savefig(fig)
pp.close()


