
import sys
from math import sqrt,sin,cos,acos,atan2,pi,asin,acos,isnan,exp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import griddata


def write_2D(fname,vec,value,error):
    with open(fname,'w') as f:
        f.write("\n".join(["\t".join(map(str,[v[1],v[2],val,err])) for v,val,err in zip(vec,value,error)]) + "\n")

def write_2D_pdf(fname,vec,value,qmax):
    mincut=len(value)/10
    vmax=np.max(value[mincut:])
    vmin=np.min(value[mincut:])
    absmin=np.min(np.abs([vmin,vmax]))
    absmax=np.max(np.abs([vmin,vmax]))
    pp = PdfPages(fname)
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)    
    xi = np.linspace(-qmax,qmax,1000)
    yi = np.linspace(-qmax,qmax,1000)
    zi = griddata((vec[:,1],vec[:,2]),value, (xi[None,:], yi[:,None]), method='cubic')

    im = plt.imshow(zi,interpolation='bilinear',origin='low',extent=[xi[0], xi[-1], yi[0], yi[-1]],norm=matplotlib.colors.SymLogNorm(absmin*0.5,vmin=-absmax,vmax=absmax))
    plt.xlabel('qy [1/nm] (laser polarization)', fontsize=16)
    plt.ylabel('qz [1/nm] (laser beam)', fontsize=16)
    print [int(-qmax+0.5)-1, int(qmax-0.5)+1, int((qmax/4)-0.5)+1]
    plt.xticks(np.arange(int(-qmax+0.5)-1, int(qmax-0.5)+1, int((qmax/4)-0.5)+1))
    plt.yticks(np.arange(int(-qmax+0.5)-1, int(qmax-0.5)+1, int((qmax/4)-0.5)+1))
    plt.tick_params(labelsize=12)
    plt.xlim(-qmax*1.05,qmax*1.05)
    plt.ylim(-qmax*1.05,qmax*1.05)
    #plt.title('Anisotropic Difference Spectra', fontsize=20)
    cbar = fig.colorbar(im,orientation='vertical',ticks=[20000,10000,5000,2000,1000,0,-1000,-2000,-5000,-10000,-20000])
    cbar.ax.tick_params(labelsize=16)
    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')   
        t.set_x(5)
    pp.savefig(fig)
    pp.close()

def write_3D(fname,vec,value):
    with open(fname,'w') as f:
        f.write("\n".join(["\t".join(map(str,[v[0],v[1],v[2],val])) for v,val in zip(vec,value)]) + "\n")

def write_xvg(fname,x_vec,y_matrix,y_error,title,legend):
    with open(fname,'w') as f:
        f.write("@    title \""+title+"\"\n")
        f.write("@    xaxis  label \"q [nm\S-1\N]\""+"\n")
        f.write("@    yaxis  label \"Intensity [e\S2\N]\""+"\n")
        f.write("@TYPE xydy"+"\n")
        for i,text in enumerate(legend):
            f.write("@    s"+str(i)+" legend  \""+ text +"\"\n")
        for y_vec,y_err in zip(y_matrix,y_error):
            f.write("@type xydy"+"\n")
            f.write("\n".join(["\t".join(map(str,[x,y,err])) for x,y,err in zip(x_vec,y_vec,y_err)]) + "\n")
            f.write("&"+"\n")

def print_mapping(outputfile,allqin,allqabsin,allqout,allqabsout,matrix,nqabs,arrows):
    outfracs=np.array([0,1.0/24,2.0/24,3.0/24,4.0/24,5.0/24,6.0/24])
    outnames=['0deg','15deg','30deg','45deg','60deg','75deg','90deg','90deg - 0deg','45deg - 0deg','45deg - 90deg','90deg - 60deg + 30deg - 0deg','0deg + 90deg - 2*45deg','azimutal (vertical)','azimutal (horizontal)','radial avg.','difference (azimutal)']
    k=nqabs
    ioutfracs=np.int_(outfracs*allqabsout.nq_out[k])
    pp = PdfPages(outputfile+"_"+str(k)+'.pdf')
    #pp_c = PdfPages(outputfile+"_curve.pdf")
    i_in = allqabsin.index[k]
    nq_in = allqabsin.nq[k]
    q = allqin[i_in:i_in+nq_in]

    theta=np.arccos(np.dot(q,arrows[0:3])/allqabsin.abs[k])
    
    plots_value=np.zeros([len(outnames),nq_in])
    plots_err=np.zeros([len(outnames),nq_in])
    
    avg=np.mean(matrix[nqabs])
    tmin=np.min(matrix[nqabs])
    tmax=np.max(matrix[nqabs])
    absmax=np.abs([tmax,-tmin])
    plt.subplots_adjust(left=0, right=1,wspace=0,top=1.15,bottom=0.15)
    q = allqin[i_in:i_in+nq_in]
    for i,iout in enumerate(ioutfracs):
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(outnames[i],y=1.1,x=0.4,fontsize = 20)
        plots_value[i]=matrix[k][iout]
        m = print_sphere(ax,q,matrix[k][iout],arrows,[-absmax,absmax],True)
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
    
    diff=matrix[k][ioutfracs[-1]]-matrix[k][ioutfracs[0]]
        
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(outnames[-9],y=1.1,x=0.4,fontsize = 20)
    plots_value[-6]=diff
    m = print_sphere(ax,q,diff,arrows,[-absmax,absmax],True)
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
    
    diff=matrix[k][ioutfracs[3]]-matrix[k][ioutfracs[0]]
        
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(outnames[-8],y=1.1,x=0.4,fontsize = 20)
    plots_value[-6]=diff
    m = print_sphere(ax,q,diff,arrows,[-absmax,absmax],True)
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

    diff=matrix[k][ioutfracs[3]]-matrix[k][ioutfracs[6]]
        
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(outnames[-7],y=1.1,x=0.4,fontsize = 20)
    plots_value[-6]=diff
    m = print_sphere(ax,q,diff,arrows,[-absmax,absmax],True)
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

    diff=matrix[k][ioutfracs[6]]-matrix[k][ioutfracs[4]]+matrix[k][ioutfracs[2]]-matrix[k][ioutfracs[0]]
        
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(outnames[-6],y=1.1,x=0.4,fontsize = 20)
    plots_value[-6]=diff
    m = print_sphere(ax,q,diff,arrows,[-absmax,absmax],True)
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

    
    diff=matrix[k][ioutfracs[-1]]+matrix[k][ioutfracs[0]]-matrix[k][ioutfracs[3]]*2
        
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(outnames[-5],y=1.1,x=0.4,fontsize = 20)
    plots_value[-5]=diff
    m = print_sphere(ax,q,diff,arrows,[-absmax,absmax],True)
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

    
    azim = integrate_matrix(matrix[k])
    diff=azim[0]-azim[1]
    absmax=np.max(np.abs(diff))
    azimmax=np.max(azim)

    for i,tazim in enumerate(azim):
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(outnames[-4+i],y=1.1,x=0.4,fontsize = 20)
        plots_value[-4+i]=tazim
        m = print_sphere(ax,q,tazim,arrows,[-azimmax,azimmax],True)
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

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(outnames[-2],y=1.1,x=0.4,fontsize = 20)
    mean=np.mean(matrix[k],axis=0)
    plots_value[-2]=mean
    m = print_sphere(ax,q,mean,arrows,[-absmax,absmax],True)
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
            
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(outnames[-1],y=1.1,x=0.4,fontsize = 20)
    plots_value[-1]=diff
    m = print_sphere(ax,q,diff,arrows,[-absmax,absmax],True)
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
    
    pp.close()

    sortindex = np.argsort(theta)    
    write_xvg(outputfile+"_"+str(k)+'.xvg',theta[sortindex],plots_value[:,sortindex],plots_err,"Weights",outnames)

def integrate_matrix(matrix):
    nq = matrix.shape[1]
    mout = np.zeros([2,nq])
    #parts = [1.0/8,3.0/8,5.0/8,7.0/8]
    #out   = [1,0,1,0]
    out = [1]*int(nq*1.0/8)+[0]*int(nq*2.0/8)+[1]*int(nq*2.0/8)+[0]*int(nq*2.0/8)+[1]*int(nq*1.0/8)
    print out
    nout = [0.0,0.0]
    print len(out)
    print matrix.shape
    for i,m in enumerate(matrix):
        mout[out[i]] += m
        nout[out[i]] += 1
    return np.array([mout[0]/nout[0],mout[1]/nout[1]])
           

def print_curve(plt,veclist,m,value):
    
    plt.plot(radius, area)


def print_sphere(ax,veclist,value,arrows,minmax,log):
    vec = np.array(veclist)
    vecT = vec.T
    print "Minmax: ",minmax
    if not minmax:
        norm=matplotlib.colors.Normalize(clip = False)
    else:
        if log:
            absmin=np.min(np.abs([minmax[0],minmax[1]]))
            absmax=np.max(np.abs([minmax[0],minmax[1]]))
        #    norm=matplotlib.colors.SymLogNorm(absmin*0.1,vmin=-absmax,vmax=absmax)
            #norm=matplotlib.colors.SymLogNorm(4,vmin=-26,vmax=26)
            norm=matplotlib.colors.Normalize(vmin=-absmax,vmax=absmax,clip = False)
            cmap=cm.coolwarm
        else:
            norm=matplotlib.colors.Normalize(vmin=minmax[0], vmax=minmax[1], clip = False)
            cmap=cm.Reds
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    m.set_array(value)
    vecnorm=np.linalg.norm(vec[0])
    vecTn=vecT/(vecnorm*1.0)
    narrows=len(arrows)/3
    colors=np.append(m.to_rgba(value),np.array([matplotlib.colors.colorConverter.to_rgba('black')]*narrows),axis=0)
    s=np.append([40]*len(value),[100]*narrows)
    ax.scatter(np.append(vecTn[0],arrows[0::3]),np.append(-vecTn[2],-arrows[2::3]),np.append(vecTn[1],arrows[1::3]),c=colors,s=s,linewidths=0.5)
    ax.set_xlim(-1.1,1.1)
    ax.set_ylim(-1.1,1.1)
    ax.set_zlim(-1.1,1.1)
    return m
