from math import sqrt,sin,cos,acos,atan2,pi,asin,acos,isnan,exp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from collections import namedtuple
from WAXSPrint import write_2D,write_2D_pdf,write_3D,write_xvg,print_mapping,print_sphere


# in eV*nm
hc = 1239.842


class ExperimentRef:
    """Scattering amplitude and intensities in the experiment reference System"""
    def __init__(self,nqabs,qabs_in,someargs):
        self.phi,self.beam,self.bulk,self.numeric,self.m,self.illustrate,self.Diff,self.polar,self.weight,self.mtype,self.legend,self.azimutal,self.scale,self.vacuum,self.nparts = someargs
       
        self.qabs = qabs_in.abs
        self.qabsin = qabs_in
        self.nqabs = nqabs
        self.nout = 0

        self.qin = []
        self.D = []
        self.l2fixed = namedtuple('l2fixed','par,hor,ver,az_hor,az_ver,avg')
        self.l2fixed.par = np.zeros([nqabs])
        self.l2fixed.hor = np.zeros([nqabs])
        self.l2fixed.ver = np.zeros([nqabs])
        self.l2fixed.az_hor = np.zeros([nqabs])
        self.l2fixed.az_ver = np.zeros([nqabs])
        self.l2fixed.avg = np.zeros([nqabs])


        self.qout = namedtuple('qout','i_out,nq_out')
        self.qout.i_out = []
        self.qout.nq_out = []
       
        add_i=self.phi
        for i in range(nqabs):
            if self.qabs[i] != 0.0:
                self.qout.i_out.append(self.qout.i_out[i-1]+add_i)
                self.qout.nq_out.append(self.phi)
                add_i=self.phi
            else:
                self.qout.i_out.append(0)
                self.qout.nq_out.append(1)
                add_i=1
            self.nout = self.nout + add_i

        self.I = np.zeros([self.nout])
        self.d_I = np.zeros([self.nout])
        self.q = np.zeros([self.nout,3])
        self.l2 = np.zeros([self.nout])
        self.D_m = np.zeros([nqabs,self.nparts])
        self.D_iso = np.zeros([nqabs,self.nparts])
         
    def init_q(self):
        q_beam = self.beam * 2. * pi * 1000 / hc
        i=0
        for j,(qabs,i_out,nq_out) in enumerate(zip(self.qabs,self.qout.i_out,self.qout.nq_out)):
            if qabs != 0:
                q_x = qabs*qabs*0.5/q_beam
                q_t = sqrt(qabs*qabs-q_x*q_x)
                fac1 =  qabs*qabs*0.25/(q_beam*q_beam)  #q^2/(4*k^2)
                print fac1
                #if not self.numeric:
                self.l2fixed.par[j] = fac1
                self.l2fixed.hor[j] = 1-fac1
                self.l2fixed.ver[j] = 0
                self.l2fixed.az_hor[j] = (pi+2)/(2*pi)*(1-fac1)
                self.l2fixed.az_ver[j] = (pi-2)/(2*pi)*(1-fac1)
                self.l2fixed.avg[j] = (1-fac1)*0.5

                for k in range(nq_out):
                    phi = k*2*pi/nq_out
                    q_y = sin(phi)*q_t
                    q_z = cos(phi)*q_t
                    q=[q_x,q_y,q_z]
                    self.q[i] = q
                    cos_ql = np.dot(q,[0,1,0])/qabs
                    self.l2[i] = cos_ql*cos_ql
                    i += 1
            else:
                self.q[i] = [0,0,0]
                self.l2[i] = 1.0/3.0
                i += 1
                self.l2fixed.par[j] = 1.0/3.0
                self.l2fixed.hor[j] = 1.0/3.0
                self.l2fixed.ver[j] = 1.0/3.0
                self.l2fixed.az_hor[j] = 1.0/3.0
                self.l2fixed.az_ver[j] = 1.0/3.0
                self.l2fixed.avg[j] = 1.0/3.0

                
    def calc_q_in_intensity_error(self,scata,scatb,scate,output):
        if not isnan(self.bulk):
            waterdensity=self.bulk
        else:
            waterdensity=0.5*(scata.solden+scatb.solden)
        self.qin = scata[0].q
        self.D = np.zeros([scata[0].nq,self.nparts])
        scale=self.scale
        scale2=scale*scale
        if self.illustrate !=0:
            pp = PdfPages(output+'_D.pdf')
        j = 0

        illuI=0
        
        for i_abs,(qabs,index,nq) in enumerate(zip(self.qabs,scata[0].qabs.index,scata[0].qabs.nq)):
            D_m,D_iso = np.zeros([self.nparts]),np.zeros([self.nparts])
            for i in range(index,index+nq):
                q = scata[0].q[i]
                self.qin.append(q)
                A = np.array([comb.sa.S[i] for comb in scata])
                A2 = np.array([comb.sa.S2[i] for comb in scata])
                if scatb == None:
                    B = np.array([0.0]*nscat)
                    B2 = np.array([0.0]*nscat)
                else:
                    B = np.array([comb.sa.S[i] for comb in scatb])
                    B2 = np.array([comb.sa.S2[i] for comb in scatb])           
                if scate != None:
                    E = np.array(scate.sa.S[i])                    
                if self.vacuum:
                    D = A2*scale2 - B2*scale2
                elif scate != None:
                    D = A2 - B2 - 2.0 * ((A-B)*E.conjugate()*waterdensity).real
                else:
                    raise NameError('Without specifing an envelope only vacuum calculation is supported.')
                self.D[j] = D
                if not self.numeric:
                    if qabs != 0.0:
                        m_q = np.dot(q,self.m[0:3])/qabs
                    else:
                        m_q = sqrt(1.0/3.0)
                    D_m = D_m + D*(3*m_q*m_q - 1)
                D_iso = D_iso + D
    
                j += 1
            
            self.D_iso[i_abs] = D_iso/nq
            self.D_m[i_abs] = D_m/nq
            
            if self.illustrate  and len(self.illustrate) > illuI and qabs >= self.illustrate[illuI] and i_abs > 0:
                illuI += 1
                tmax=np.max(np.mean(self.D[index:index+nq],axis=1))/1000
                tmin=np.min(np.mean(self.D[index:index+nq],axis=1))/1000
                fig = plt.figure(figsize=(7, 6))
                ax = fig.add_subplot(111, projection='3d')
                ax.set_title("|q| = %8.1f"%qabs,y=1.05,x=0.4,fontsize = 20)
                m=print_sphere(ax,scata[0].q[index:index+nq],np.mean(self.D[index:index+nq],axis=1)/1000,np.array(self.m),[-tmin,tmax],True)
                for a in ax.w_xaxis.get_ticklabels()+ax.w_yaxis.get_ticklabels()+ax.w_zaxis.get_ticklabels():
                    a.set_visible(False) 
                cax = fig.add_axes([0.05,0.3,0.2,0.4])
                cax.get_xaxis().set_visible(False)
                cax.get_yaxis().set_visible(False)
                cax.patch.set_alpha(0)
                cax.set_frame_on(False)
                cbar = fig.colorbar(m,orientation='vertical',ax=cax,ticks=plt.MaxNLocator(nbins=6))
                cbar.ax.tick_params(labelsize=20) 
                pp.savefig(fig)
        if self.illustrate !=0:
            pp.close()

    def combine_intensities(self,intensities,populations):
        self.D = np.zeros(intensities[0].D.shape)
        for intensity,population in zip(intensities,populations):
            self.D_iso += intensity.D_iso * population
            self.D_m += intensity.D_m * population
            self.D += intensity.D * population      

    def do_passiv_rot(self,trot):
        if self.polar == 'circular':
            print 'Two print 2D patterns with circular polarisation is currently not supported.'
            return
        if self.weight != "cos2":
            raise NameError('Weight not supported. Use numeric version.')
        difffactor=exp(-(trot*self.Diff*6))
        for i_abs,(i_out,nq_out,D_m,D_iso) in enumerate(zip(self.qout.i_out,self.qout.nq_out,self.D_m,self.D_iso)):
            for i_q_exp in range(i_out,i_out+nq_out):
                l2 = self.l2[i_q_exp]
                if self.polar == "linear":
                    l2 = self.l2[i_q_exp]
                elif self.polar == "circular":
                    l2 = (self.l2[i_q_exp] + self.l2fixed.par[i_abs])/2
                temp = (l2*1.5-0.5)*difffactor
                temp2 = np.array([temp]*self.nparts).T
                temp3 = temp2*D_m
                if self.mtype == 'linear':
                    I = (temp3+D_iso)
                elif self.mtype == 'planar':
                    I = (-temp3*0.5+D_iso)
                self.I[i_q_exp] = np.mean(I)
                self.d_I[i_q_exp] = np.std(I)/sqrt(self.nparts)

             
    def do_passiv_rot_numeric(self,function,trot):
        print "Do numeric mapping."
        if self.polar == 'circular':
            print 'Two print 2D patterns with circular polarisation is currently not supported.'
            return
        difffactor=exp(-(trot*self.Diff*6))
        transmatrix = function(self.qin,self.qabsin,self.q,self.qout)
        print "Transition matrix stored."
        for k,(i_in,nq_in) in enumerate(zip(self.qabsin.index,self.qabsin.nq)):
            if nq_in > 1:
                temp = np.tensordot(self.D[i_in:i_in+nq_in],transmatrix[k],axes=([0],[1])).T/nq_in
                Itemp=np.append(Itemp,temp,axis=0)
            else:
                if k != 0:
                    raise NameError("This should not happen. (do_passiv_rot_numeric)")
                Itemp=np.array([self.D[k]])
        I = Itemp
        print "Transition done."
        self.I = np.mean(I,axis=1)
        self.d_I = np.std(I,axis=1)/sqrt(self.nparts)
                

    def I_of_lq2(self,l2,difffusion,trot):
        if difffusion:
            difffactor=exp(-(trot*self.Diff*6))
        else:
            difffactor=1
        if self.polar == "circular":
            l2 = (l2 + self.l2fixed.par)/2
        temp = (l2*1.5-0.5)*difffactor
        temp = np.array([temp]*self.nparts).T
        temp = temp*self.D_m
        if self.mtype == 'linear':
            I = (temp+self.D_iso)
        elif self.mtype == 'planar':
            I = (-temp*0.5+self.D_iso)
        return (np.mean(I,axis=1),np.std(I,axis=1)/sqrt(self.nparts))
        
                        
    def write_intensity(self,trot,output):
        if not self.polar == 'circular':
            write_2D(output+"_aniso.dat",self.q,self.I,self.d_I)
            write_2D_pdf(output+'_aniso.pdf',self.q,self.I,self.qabs[-1])
        write_xvg(output+"_iso.xvg",self.qabs,[np.mean(self.D_iso,axis=1)],[np.std(self.D_iso,axis=1)/sqrt(self.nparts)],"Difference WAXS pattern: Isotrop",[self.legend])
        write_xvg(output+"_alliso.xvg",self.qabs,self.D_iso.T,np.zeros(self.D_iso.shape).T,"WAXS pattern: Isotrop",map(str,range(0,self.D_iso.shape[1])))
        write_xvg(output+"_m.xvg",self.qabs,[np.mean(self.D_m,axis=1)],[np.std(self.D_m,axis=1)/sqrt(self.nparts)],"Difference WAXS pattern: Anisotrop",[self.legend])

        parallel,Errparallel=self.I_of_lq2(self.l2fixed.par,True,trot)
        write_xvg(output+"_par.xvg",self.qabs,[parallel],[Errparallel],"Difference WAXS pattern: Parallel laser polarisation.",[self.legend])
        
        if self.azimutal:
            Ihorizontal,Errhorizontal=self.I_of_lq2(self.l2fixed.az_hor,True,trot)
            Ivertical,Errvertical=self.I_of_lq2(self.l2fixed.az_ver,True,trot)
            temp = np.array([self.l2fixed.az_hor-self.l2fixed.az_ver]*self.nparts).T
        else:
            Ihorizontal,Errhorizontal=self.I_of_lq2(self.l2fixed.hor,True,trot)
            Ivertical,Errvertical=self.I_of_lq2(self.l2fixed.ver,True,trot)
            temp = np.array([self.l2fixed.hor-self.l2fixed.ver]*self.nparts).T
        if self.mtype == 'linear':
            HminV = temp*self.D_m*exp(-(trot*self.Diff*6))*1.5
        elif self.mtype == 'planar':
            HminV = -temp*self.D_m*exp(-(trot*self.Diff*6))*1.5/2
        if self.polar == "circular":
            HminV = HminV/2


        write_xvg(output+"_V.xvg",self.qabs,[Ivertical],[Errvertical],"WAXS pattern: Vertical",[self.legend])
        write_xvg(output+"_H.xvg",self.qabs,[Ihorizontal],[Errhorizontal],"WAXS pattern: Horizontal",[self.legend])
        write_xvg(output+"_HminV.xvg",self.qabs,[np.mean(HminV,axis=1)],[np.std(HminV,axis=1)/sqrt(self.nparts)],"WAXS pattern: Horizontal - Vertical",[self.legend])
        
        Iavg,Erravg=self.I_of_lq2(self.l2fixed.avg,True,trot)
        write_xvg(output+"_avg.xvg",self.qabs,[Iavg],[Erravg],"WAXS pattern: Rotational Averaged",[self.legend])


    def add_plots(self,ax_avg,ax_HminV,num,trot,shift_avg,shift_HminV,label):
        if self.azimutal:
            Ihorizontal,Errhorizontal=self.I_of_lq2(self.l2fixed.az_hor,True,trot)
            Ivertical,Errvertical=self.I_of_lq2(self.l2fixed.az_ver,True,trot)
            temp = np.array([self.l2fixed.az_hor-self.l2fixed.az_ver]*self.nparts).T
        else:
            Ihorizontal,Errhorizontal=self.I_of_lq2(self.l2fixed.hor,True,trot)
            Ivertical,Errvertical=self.I_of_lq2(self.l2fixed.ver,True,trot)
            temp = np.array([self.l2fixed.hor-self.l2fixed.ver]*self.nparts).T
        if self.mtype == 'linear':
            HminV = temp*self.D_m*exp(-(trot*self.Diff*6))*1.5
        elif self.mtype == 'planar':
            HminV = -temp*self.D_m*exp(-(trot*self.Diff*6))*1.5/2
        if self.polar == "circular":
            HminV = HminV/2
        Iavg,Erravg=self.I_of_lq2(self.l2fixed.avg,True,trot)

        temp = self.qabs*Iavg
        if shift_avg == None:
            shift_avg = np.max(temp)*0.5
            line_avg = [shift_avg]*len(self.qabs)
            ax_avg.plot(self.qabs,line_avg,'--k')
        y=temp-shift_avg*num
        yerr=self.qabs*Erravg
        line_avg = [-shift_avg*num]*len(self.qabs)
        ax_avg.plot(self.qabs,line_avg,'--k')
        ax_avg.plot(self.qabs,y,'k')
        ax_avg.fill_between(self.qabs, y-yerr, y+yerr,facecolor='grey',edgecolor="none")
        ax_avg.text((self.qabs[-1]-self.qabs[0])*0.8+self.qabs[0],shift_avg*(-num+0.4),label)

        temp = self.qabs*np.mean(HminV,axis=1)
        if shift_HminV == None:
            shift_HminV = np.max([np.absolute(np.min(temp)),np.max(temp)])
            line_HminV = [shift_HminV]*len(self.qabs)
            ax_HminV.plot(self.qabs,line_HminV,'--k')
        y=temp-shift_HminV*num
        yerr=self.qabs*np.std(HminV,axis=1)/sqrt(self.nparts)
        line_HminV = [-shift_HminV*num]*len(self.qabs)
        ax_HminV.plot(self.qabs,line_HminV,'--k')
        ax_HminV.plot(self.qabs,y,'k')
        ax_HminV.fill_between(self.qabs, y-yerr, y+yerr,facecolor='grey',edgecolor="none")
        ax_HminV.text((self.qabs[-1]-self.qabs[0])*0.8+self.qabs[0],shift_HminV*(-num+0.4),label)

        return shift_avg,shift_HminV


    def integrate(self,start,end,mirrow):
        Iout = np.zeros([self.nqabs,self.nparts])
        for k,(nq,i_out) in enumerate(zip(self.qout.nq_out,self.qout.i_out)):
            if nq == 1:
                Iout[k] = self.I[i_out]
            else:
                lim=np.int_(np.array([start,end])*nq)
                print lim
                Iout[k] = np.mean(self.I[i_out+lim[0]:i_out+lim[1]],axis=0)
                if mirrow:
                    half=int(nq/2)
                    lim2=np.int_([lim[0]+half,lim[1]+half]*nq)
                    lim2=np.mod(lim2,nq)
                    Iout[k] = (Iout[k] + np.mean(self.I[i_out+lim2[0]:i_out+lim2[1]],axis=0))*0.5
        return Iout

    def pick(self,start,mirrow):
        Iout = np.zeros([self.nqabs,self.nparts])
        for k,(nq,i_out) in enumerate(zip(self.qout.nq_out,self.qout.i_out)):
            if nq == 1:
                Iout[k] = self.I[i_out]
            else:
                lim=np.int_(start*nq)
                Iout[k] = self.I[i_out+lim]
                if mirrow:
                    half=int(nq/2)
                    lim2=np.int_(lim+half)
                    lim2=np.mod(lim2,nq)
                    Iout[k] = (Iout[k] + self.I[i_out+lim2])*0.5
        return Iout
                
    def write_numeric(self,output):
        write_2D(output+"_aniso.dat",self.q,self.I,self.d_I)
        write_2D_pdf(output+'_aniso.pdf',self.q,self.I,self.qabs[-1])
        write_xvg(output+"_iso.xvg",self.qabs,[np.mean(self.D_iso,axis=1)],[np.std(self.D_iso,axis=1)/sqrt(self.nparts)],"Difference WAXS pattern: Isotrop",[self.legend])
        if self.azimutal:
            allIhorizontal=self.integrate(0.125,0.375,False)
            allIvertical=self.integrate(0.375,0.625,False)
            allI60=self.integrate(1.0/6-1.0/8,1.0/6+1.0/8,False)
            allI45=self.integrate(1.0/8-1.0/8,1.0/8+1.0/8,False)
            allI30=self.integrate(1.0/12-1.0/8,1.0/12+1.0/8,False)
        else:
            allIhorizontal=self.pick(0.25,True)
            allIvertical=self.pick(0,True)
            allI60=self.pick(1.0/6,True)
            allI45=self.pick(1.0/8,True)
            allI30=self.pick(1.0/12,True)
        
        Ihorizontal,Errhorizontal=np.mean(allIhorizontal,axis=1),np.std(allIhorizontal,axis=1)/sqrt(self.nparts)
        Ivertical,Errvertical=np.mean(allIvertical,axis=1),np.std(allIvertical,axis=1)/sqrt(self.nparts)
        diff=allIhorizontal-allIvertical
        IHminV,ErrHminV=np.mean(diff,axis=1),np.std(diff,axis=1)/sqrt(self.nparts)

        diff45_0=allI45-allIvertical
        I45_0,Err45_0=np.mean(diff45_0,axis=1),np.std(diff45_0,axis=1)/sqrt(self.nparts)
        diff45_90=allI45-allIhorizontal
        I45_90,Err45_90=np.mean(diff45_90,axis=1),np.std(diff45_90,axis=1)/sqrt(self.nparts)
        diff90_45_0=allIhorizontal+allIvertical-allI45*2
        I90_45_0,Err90_45_0=np.mean(diff90_45_0,axis=1),np.std(diff90_45_0,axis=1)/sqrt(self.nparts)
        diff90_60_30_0=allIhorizontal-allI60+allI30-allIvertical
        I90_60_30_0,Err90_60_30_0=np.mean(diff90_60_30_0,axis=1),np.std(diff90_60_30_0,axis=1)/sqrt(self.nparts)
        

        write_xvg(output+"_V.xvg",self.qabs,[Ivertical],[Errvertical],"WAXS pattern: Vertical",[self.legend])
        write_xvg(output+"_H.xvg",self.qabs,[Ihorizontal],[Errhorizontal],"WAXS pattern: Horizontal",[self.legend])
        write_xvg(output+"_HminV.xvg",self.qabs,[IHminV],[ErrHminV],"WAXS pattern: Horizontal - Vertical",[self.legend])
        write_xvg(output+"_45_0.xvg",self.qabs,[I45_0],[Err45_0],"WAXS pattern: 45 deg - 0 deg",[self.legend])
        write_xvg(output+"_45_90.xvg",self.qabs,[I45_90],[Err45_90],"WAXS pattern: 45 deg - 90 deg",[self.legend])
        write_xvg(output+"_90_45_0.xvg",self.qabs,[I90_45_0],[Err90_45_0],"WAXS pattern: 0 deg + 90 deg - 45 deg",[self.legend])
        write_xvg(output+"_90_60_30_0.xvg",self.qabs,[I90_60_30_0],[Err90_60_30_0],"WAXS pattern: 90 deg - 60deg + 30 deg - 0 deg",[self.legend])
        allIavg=self.integrate(0,1,False)
        Iavg,Erravg=np.mean(allIavg,axis=1),np.std(allIavg,axis=1)/sqrt(self.nparts)
        write_xvg(output+"_avg.xvg",self.qabs,[Iavg],[Erravg],"WAXS pattern: Radial Averaged",[self.legend])      
