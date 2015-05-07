#! /usr/bin/env python

import sys
import argparse
from argparse import RawTextHelpFormatter
from math import sqrt,sin,cos,acos,atan2,pi,asin,acos,isnan,exp
import numpy as np
import scipy 


from mypython.WAXSPrint import print_mapping
from mypython.ExperimentRef import ExperimentRef
from mypython.ProteinRef import ProteinRef,check_consistency

# in eV*nm
hc = 1239.842

            
def doublecos2weight(qin,qout):
    raise NameError('doublecos2weight: Still to be done.')


def integrate_cos2(a,lq,mq):
    psi = np.arange(0,2*np.pi,np.pi/50)
    cs = np.cos(psi)
    lt=np.sqrt(1-lq*lq)
    mt=np.sqrt(1-mq*mq)
    f =  1-np.exp(-a*(lq*lq*mq*mq+2*lq*mq*lt*mt*cs+lt*lt*mt*mt*cs*cs))
    return np.mean(f)  

def addcos2_matrix(allqin,allqabsin,allqout,allqabsout):
    difffactor=exp(-(args.t*args.D*6))
    output = np.zeros([len(allqabsin.abs),args.phi,allqabsin.nq[-1]])
    illuI = 0
    int1=np.abs(args.w_args[0])
    int2=np.abs(args.w_args[2])
    totex1=1-(1*np.sqrt(np.pi)*scipy.special.erf(np.sqrt(int1)))/(2*np.sqrt(int1)) if int1 != 0 else 0
    totex2=1-(1*np.sqrt(np.pi)*scipy.special.erf(np.sqrt(int2)))/(2*np.sqrt(int2)) if int2 != 0 else 0
    norm1=1.0/totex1
    norm2=1.0/totex2
    refint1=args.w_args[1]
    refint2=args.w_args[3]
    print totex1,totex2
    for k,(qabs,i_in,nq_in,i_out,nq_out) in enumerate(zip(allqabsin.abs,allqabsin.index,allqabsin.nq,allqabsout.i_out,allqabsout.nq_out)):
        print "Calculation for qabs = "+str(qabs)
        if qabs != 0:
            for i,q_out in enumerate(allqout[i_out:i_out+nq_out]):
                lq=np.dot(q_out,[0,1,0])/qabs
                for j,q_in in enumerate(allqin[i_in:i_in+nq_in]):
                    mq=np.dot(q_in,args.m[0:3])/qabs
                    if args.mtype == 'linear':
                        output[k][i][j] = refint1*norm1*integrate_cos2(int1,lq,mq)+refint2*norm2*integrate_cos2(int2,lq,mq)
                    elif args.mtype == 'planar':
                        raise NameError('planar: Still to be done.')
        else:
            output[k][0][0] = refint1+refint2
     
        if args.illustrate and len(args.illustrate) > illuI and qabs >= args.illustrate[illuI] and k > 0:
            illuI += 1
            print_mapping(args.o+'_l2_'+str(qabs),allqin,allqabsin,allqout,allqabsout,output,k,np.array(args.m))

    return output


def cos2_matrix(allqin,allqabsin,allqout,allqabsout):
    difffactor=exp(-(args.t*args.D*6))
    output = np.zeros([len(allqabsin.abs),args.phi,allqabsin.nq[-1]])
    illuI = 0
    for k,(qabs,i_in,nq_in,i_out,nq_out) in enumerate(zip(allqabsin.abs,allqabsin.index,allqabsin.nq,allqabsout.i_out,allqabsout.nq_out)):
        print "Calculation for qabs = "+str(qabs)
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
       
        if args.illustrate and len(args.illustrate) > illuI and qabs >= args.illustrate[illuI] and k > 0:
            illuI += 1
            print_mapping(args.o+'_l2_'+str(qabs),allqin,allqabsin,allqout,allqabsout,output,k,np.array(args.m))
    return output

def dcos2_matrix(allqin,allqabsin,allqout,allqabsout):
    output = np.zeros([len(allqabsin.abs),args.phi,allqabsin.nq[-1]])
    m=np.array(args.m[0:3])
    n=np.array(args.m[3:6])
    if args.mtype == 'planar':
        raise NameError('Planar transition moment not supported for two transition moments.')
    for k,(qabs,i_in,nq_in,i_out,nq_out) in enumerate(zip(allqabsin.abs,allqabsin.index,allqabsin.nq,allqabsout.i_out,allqabsout.nq_out)):
        print "Calculation for qabs = "+str(qabs)
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
        if args.illustrate !=0 and k == 1:
            print_mapping(args.o+'_l2',allqin,allqabsin,allqout,allqabsout,output,k,np.array(args.m))

    return output

funcdict = {
  'cos2': cos2_matrix,
  'dcos2': dcos2_matrix,
  'addcos2': addcos2_matrix
}
    
parser = argparse.ArgumentParser(description='Anisotropic pattern by passive rotation.\n\n'
                                 'Input files:\n'
                                 '*_averagA.xvg for system excited and in ground state.\n\n'
                                 'Output files:\n'
                                 '_aniso.dat / _aniso.pdf       2D Scattering pattern.\n'
                                 '_iso.xvg                      Isotropic component.\n'
                                 '_m.xvg                        Anisotropic component.\n'
                                 '_avg.xvg                      Azimutal average.\n'
                                 '_par.xvg                      Laser polarisation parallel.\n'
                                 '_H.xvg                        Horizontal cut / 90deg average.\n'
                                 '_V.xvg                        Vertical cut / 90deg average.\n'
                                 '_HminV.xvg                    Difference between horizontal and vertical.\n'
                                 '_alliso.xvg                   Isotropic component from all bins individually.\n'
                                 '_D.pdf                        Difference term for selected momenta.\n\n'
                                 'Additional output for numeric calculations:\n'
                                 '_l2_###.xvg                   Weigts for caculation of different quantities.\n'
                                 '_45_0.xvg / _45_90.xvg / _90_45_0.xvg / _90_60_30_0.xvg   additional quantities.\n'
                                 
                                 , formatter_class=RawTextHelpFormatter)
parser.add_argument('-fa', nargs='+', metavar='_averagA.xvg',help='Average scattering amplitude of System A.')
parser.add_argument('-fb', nargs='+', metavar='_averagA.xvg',help='Average scattering amplitude of System B.')
parser.add_argument('-env', metavar='',help='Fourier transform of the envelope.')
parser.add_argument('-env_ab', action='store_true', help='Assume envelope files for input -fa and -fb.')
parser.add_argument('-o', metavar='.xvg',help='Output.')
parser.add_argument('-vacuum', action='store_true', help='Do vacuum calculation.')
parser.add_argument('-phi', type=int, default=360, help='Number of azimutal data points in 2D output.')
parser.add_argument('-weight', choices=['cos2','dcos2','addcos2'], default='cos2', help='Weight of different orientations.')
parser.add_argument('-m', type=float, default=[1,0,0],metavar="", nargs='+', help='Vector of excitation moment(s). (cos2,addcos2: one vector; dcos2: two vector)')
parser.add_argument('-w_args', type=float, nargs=4,metavar=('intensity1','weight1','intensity2','weight2'), help='(only used for addcos2) laser intensity and corresponding weight after normalisation.\n')
parser.add_argument('-mtype', choices=['linear','planar'], default='linear',
                   help='Type of the excitation moment.')
parser.add_argument('-polar', choices=['linear','circular'], default='linear',
                   help='Type of excitation laser polarisation.')
parser.add_argument('-beam', type=float, default=12.0, help='X-ray beam energy. [keV]')
parser.add_argument('-bulk', type=float, default='nan', help='Set bulk water density. [e/nm^-3]')
parser.add_argument('-D', type=float, default=0.0, help='Rotational diffusion constant. [1/(unit of time)]')
parser.add_argument('-t', type=float, default=0.0, help='Delay time. [unit of time]')
parser.add_argument('-numeric', action='store_true',help='Do numeric integration over SO(3) (slow), do not use two components. Numeric is required for weight=dcos2/addcos2.')
parser.add_argument('-legend', help='Name in legend for plot.')
parser.add_argument('-azimutal', action='store_true', help='Calculate horizontal and vertical 90deg azimutal averages instead of cuts.')
parser.add_argument('-scale', type=float, default=1, help="Scaling of the electron density.")
parser.add_argument('-illustrate', nargs='+', type=float, help="Illustrate the rotational averaging. Enter q values to illustrate.")


args = parser.parse_args()


nparts=len(args.fa)

if args.vacuum:
    raise NameError('Vacuum calculations are not approved. Check code first.')

if args.polar == 'circular':
    raise NameError('Calculations with circular exciation laser polarisation are not approved. Check code first.')

allA = [ProteinRef() for fa in args.fa]
for fname,A in zip(args.fa,allA):
    if not args.env_ab:
        A.read_file(fname)
    else:
        A.read_envelope(fname)
    check_consistency(allA[0],A)

allB = [ProteinRef() for fb in args.fb]
if len(args.fa) != len(args.fb):
    raise NameError('The same number of files has to be given as an input for A and B.')
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
I=ExperimentRef(allA[0].nqabs,allA[0].qabs,someargs)
if args.env:
    I.calc_q_in_intensity_error(allA,allB,E,args.o)
else:
    I.calc_q_in_intensity_error(allA,allB,None,args.o)
I.init_q()
if not args.numeric:
    I.do_passiv_rot(args.t)
    I.write_intensity(args.t,args.o)
else:
    I.do_passiv_rot_numeric(funcdict[args.weight],args.t)
    I.write_numeric(args.o)

