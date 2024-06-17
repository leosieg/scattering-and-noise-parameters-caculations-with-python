#
import numpy as np
import matplotlib.pyplot as plt
import lib_func.lib_smith_stab as sms
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def smith_z(val):
    #
    # Smith chart impedance(z) version
    #
    # val -> values of real and imaginary impedance part circles
    #        example: np.array([0., 0.2, 0.5, 1., 2., 5.])
    #                     
    pk = 200
    phi = np.linspace(0,2*np.pi,pk)
    #
    counts = val
    L_val = counts.size
    #
    q = 0
    # circles constant real part
    while q <= L_val-1:
        #
        x0 = val[q]/(val[q]+1)
        ra = 1./(val[q]+1)
        #
        z = x0 + ra* np.exp(phi*1.j)
        #
        if ra == 1. or ra == 0.5:
            plt.plot(np.real(z),np.imag(z),color='r',
                     lw=2)
        else:
            plt.plot(np.real(z),np.imag(z),color='r',
                     lw=1)
        q = q+1 
    #
    # match circle z -> y g = 1
    z = -0.5 + 0.5*np.exp(phi*1.j);
    plt.plot(np.real(z),np.imag(z),color='m',
             linewidth=1)
    #
    x0 = 1;
    q = 1;
    # circles constant -imaginary part
    while q <= L_val-1:
        #
        y0 = 1./-val[q]
        ra = 1./np.abs(val[q])
        #
        fis = np.pi/2.
        fie = np.pi-np.arccos(2.*val[q]/(1+val[q]**2))
        #
        if np.abs(val[q]) > 1:
            #
            fie = np.pi+np.arccos(2.*val[q]/(1+val[q]**2))
            #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j);
        #
        plt.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    q = 1
    # circles constant +imaginary part
    while q <= L_val-1:
        #
        y0 = 1./val[q]
        ra = 1./np.abs(val[q])
        #
        fis = -np.arccos(-2.*val[q]/(1+val[q]**2))
        if np.abs(val[q]) > 1:
            #
            fis = -2.*np.pi+np.arccos(-2.*val[q]/(1+val[q]**2));
            #
        fie = -np.pi/2.
        #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j)
        #
        plt.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    # real axis
    xi = np.arange(-1, 1, 0.1)
    yi = 0*xi
    plt.plot(xi,yi,color='b',lw=2)
    #
    q = 0
    # text real axis
    while q <= L_val-1:
        xq = (val[q]-1)/(val[q]+1)+0.02;
        plt.text(xq,0.05,str(val[q]),color='r',
                 fontweight='bold',fontsize=8)
        q = q+1;
    #
    q = 0
    # text +imag
    while q <= L_val-1:
        #
        xq = (1.j*val[q]-1.)/(1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '+j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            plt.text(xp,yp,sjn,color='b',va='bottom',ha='center',
                         fontweight='bold',fontsize=8)
        elif xp < 0:
            #
            plt.text(xp,yp,sjn,color='b',va='bottom',ha='right',
                         fontweight='bold',fontsize=8)
        else:
            plt.text(xp,yp,sjn,color='b',va='bottom',ha='left',
                         fontweight='bold',fontsize=8)
        q = q+1;
    #
    q = 1
    # text -imag
    while q <= L_val-1:
        #
        xq = (-1.j*val[q]-1.)/(-1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '-j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            plt.text(xp,yp,sjn,color='b',va='top',ha='center',
                         fontweight='bold',fontsize=8)
        elif xp < 0:
            #
            plt.text(xp,yp,sjn,color='b',va='top',ha='right',
                         fontweight='bold',fontsize=8)
        else:
            plt.text(xp,yp,sjn,color='b',va='top',ha='left',
                         fontweight='bold',fontsize=8)
        q = q+1;
    #
#
#####################################################################    
#
def smith_y(val):
    #
    # Smith chart admittance(y) version
    #
    # val -> values of real and imaginary admittance part circles
    #        example: np.array([0., 0.2, 0.5, 1., 2., 5.])
    #   
    pk = 200
    phi = np.linspace(0,2*np.pi,pk)
    #
    counts = val
    L_val = counts.size
    #
    q = 0
    # circles constant real part
    while q <= L_val-1:
        #
        x0 = -val[q]/(val[q]+1)
        ra = 1./(val[q]+1)
        #
        z = x0 + ra* np.exp(phi*1.j)
        #
        if ra == 1. or ra == 0.5:
            plt.plot(np.real(z),np.imag(z),color='r',
                     lw=2)
        else:
            plt.plot(np.real(z),np.imag(z),color='r',
                     lw=1)
        q = q+1 
    #
    # match circle z -> y g = 1
    z = 0.5 + 0.5*np.exp(phi*1.j);
    plt.plot(np.real(z),np.imag(z),color='m',
             linewidth=1)
    #
    x0 = -1;
    q = 1;
    # circles constant -imaginary part
    while q <= L_val-1:
        #
        y0 = 1./-val[q]
        ra = 1./np.abs(val[q])
        #
        fis = np.pi/2.
        fie = np.pi-np.arccos(-2.*val[q]/(1+val[q]**2))
        #
        if np.abs(val[q]) > 1:
            #
            fie = -np.pi+np.arccos(-2.*val[q]/(1+val[q]**2))
            #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j);
        #
        plt.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    q = 1
    # circles constant +imaginary part
    while q <= L_val-1:
        #
        y0 = 1./val[q]
        ra = 1./np.abs(val[q])
        #
        fis = -np.arccos(2.*val[q]/(1+val[q]**2))
        if np.abs(val[q]) > 1:
            #
            fis = np.pi-np.arccos(-2.*val[q]/(1+val[q]**2));
            #
        fie = -np.pi/2.
        #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j)
        #
        plt.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    # real axis
    xi = np.arange(-1, 1.1, 0.1)
    yi = 0*xi
    plt.plot(xi,yi,color='b',lw=2)
    #
    q = 0
    # text real axis
    while q <= L_val-1:
        xq = -(val[q]-1)/(val[q]+1)+0.02;
        plt.text(xq,0.05,str(val[q]),color='r',
                 fontweight='bold',fontsize=8)
        q = q+1;
    #
    q = 0
    # text +imag
    while q <= L_val-1:
        #
        xq = -(1.j*val[q]-1.)/(1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '+j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            plt.text(xp,yp,sjn,color='b',va='top',ha='center',
                         fontweight='bold',fontsize=8)
        elif xp < 0:
            #
            plt.text(xp,yp,sjn,color='b',va='top',ha='right',
                         fontweight='bold',fontsize=8)
        else:
            plt.text(xp,yp,sjn,color='b',va='top',ha='left',
                         fontweight='bold',fontsize=8)
        q = q+1;
    #
    q = 1
    # text -imag
    while q <= L_val-1:
        #
        xq = -(-1.j*val[q]-1.)/(-1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '-j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            plt.text(xp,yp,sjn,color='b',va='bottom',ha='center',
                         fontweight='bold',fontsize=8)
        elif xp < 0:
            #
            plt.text(xp,yp,sjn,color='b',va='bottom',ha='right',
                         fontweight='bold',fontsize=8)
        else:
            plt.text(xp,yp,sjn,color='b',va='bottom',ha='left',
                         fontweight='bold',fontsize=8)
        q = q+1;
    #
#
#####################################################################
#    
def ssks(S,dsp):
    #
    # stability circle source plane
    #
    #   S -> S-matrix
    # dsp -> display paramter for Smith chart
    #         'z' -> impedance version, 'y' -> admittance version
    # 
    s11 = S[0,0]
    s12 = S[0,1]
    s21 = S[1,0]
    s22 = S[1,1]
    #
    ds = s11*s22-s12*s21
    c1 = s11-np.conj(s22)*ds
    c2 = s22-np.conj(s11)*ds
    #
    # stability factors
    mu1 = (1.-np.abs(s11)**2)/(np.abs(c2)+np.abs(s12*s21));
    mu2 = (1.-np.abs(s22)**2)/(np.abs(c1)+np.abs(s12*s21));
    #
    mu = np.array([mu1, mu2])
    #
    # stability circle
    cs = np.conj(c1)/(np.abs(s11)**2-np.abs(ds)**2)
    #
    rs = np.abs(s12*s21/(np.abs(s11)**2-np.abs(ds)**2))
    #
    # point density
    pk = 200
    phi = np.linspace(0.,2.*np.pi,pk)
    #
    z = cs + rs*np.exp(1.j*phi);
    #
    val_smith = np.array([0., 0.2, 0.5, 1., 2., 5.])
    #
    # Smith chart
    if dsp == 'z':
       sms.smith_z(val_smith)
    else:
       sms.smith_y(val_smith) 
    #
    plt.axis('off')
    plt.axis('equal')
    #
    plt.plot(np.real(z), np.imag(z), color = 'k', lw = 2)
    #
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    #
    # s11
    xcs = np.real(s11);
    ycs = np.imag(s11);
    #
    plt.plot(xcs, ycs, marker=".", markersize=5)
    #
    plt.text(xcs,ycs,'$s_{11}$',color='k',va='bottom',ha='center',
                         fontweight='bold',fontsize=10)
    #
    plt.text(-1,1.1,'stab. circle, source plane',color='k',
             fontweight='bold',fontsize=10)
    #
    print(' ')
    print('stab. factors mu = ', mu.round(2))
    #
#
#####################################################################
#    
def sskl(S,dsp):
    #
    # stability circle load plane 
    #
    #   S -> S-matrix
    # dsp -> display paramter for Smith chart
    #         'z' -> impedance version, 'y' -> admittance version
    #
    s11 = S[0,0]
    s12 = S[0,1]
    s21 = S[1,0]
    s22 = S[1,1]
    #
    ds = s11*s22-s12*s21
    c1 = s11-np.conj(s22)*ds
    c2 = s22-np.conj(s11)*ds
    #
    # stability factors
    mu1 = (1.-np.abs(s11)**2)/(np.abs(c2)+np.abs(s12*s21));
    mu2 = (1.-np.abs(s22)**2)/(np.abs(c1)+np.abs(s12*s21));
    #
    mu = np.array([mu1, mu2])
    #
    # stability circle
    cl = np.conj(c2)/(np.abs(s22)**2-np.abs(ds)**2)
    #
    rl = np.abs(s12*s21/(np.abs(s22)**2-np.abs(ds)**2))
    #
    # point density
    pk = 200
    phi = np.linspace(0.,2.*np.pi,pk)
    #
    z = cl + rl*np.exp(1.j*phi);
    #
    val_smith = np.array([0., 0.2, 0.5, 1., 2., 5.])
    #
    # Smith chart
    if dsp == 'z':
       sms.smith_z(val_smith)
    else:
       sms.smith_y(val_smith) 
    #
    plt.axis('off')
    plt.axis('equal')
    #
    plt.plot(np.real(z), np.imag(z), color = 'k', lw = 2.)
    #
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    #
    # s11
    xcs = np.real(s22);
    ycs = np.imag(s22);
    #
    plt.plot(xcs, ycs, marker=".", markersize=5)
    #
    plt.text(xcs,ycs,'$s_{22}$', color='k', va='top', ha='center',
                         fontweight='bold', fontsize = 10)
    #
    plt.text(-1,1.1,'stab. circle, load plane',color = 'k',
             fontweight = 'bold', fontsize = 10)
    #
    print('stab. factors mu = ', mu.round(2))
    #
#
#####################################################################
# 
def smith3_z(val,ax):
    #
    # Smith chart impedance(z) version for 3d projection
    #
    # val -> values of real and imaginary impedance part circles
    #        example: np.array([0., 0.2, 0.5, 1., 2., 5.])
    #  ax -> axis variable of 3d projection 
    pk = 200
    phi = np.linspace(0,2*np.pi,pk)
    #
    counts = val
    L_val = counts.size
    #
    q = 0
    # circles constant real part
    while q <= L_val-1:
        #
        x0 = val[q]/(val[q]+1)
        ra = 1./(val[q]+1)
        #
        z = x0 + ra* np.exp(phi*1.j)
        #
        if ra == 1. or ra == 0.5:
            ax.plot(np.real(z),np.imag(z),color='r',
                     lw=2)
        else:
            ax.plot(np.real(z),np.imag(z),color='r',
                     lw=1)
        q = q+1 
    #
    # match circle z -> y g = 1
    z = -0.5 + 0.5*np.exp(phi*1.j);
    ax.plot(np.real(z),np.imag(z),color='m',
             linewidth=1)
    #
    x0 = 1;
    q = 1;
    # circles constant -imaginary part
    while q <= L_val-1:
        #
        y0 = 1./-val[q]
        ra = 1./np.abs(val[q])
        #
        fis = np.pi/2.
        fie = np.pi-np.arccos(2.*val[q]/(1+val[q]**2))
        #
        if np.abs(val[q]) > 1:
            #
            fie = np.pi+np.arccos(2.*val[q]/(1+val[q]**2))
            #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j);
        #
        plt.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    q = 1
    # circles constant +imaginary part
    while q <= L_val-1:
        #
        y0 = 1./val[q]
        ra = 1./np.abs(val[q])
        #
        fis = -np.arccos(-2.*val[q]/(1+val[q]**2))
        if np.abs(val[q]) > 1:
            #
            fis = -2.*np.pi+np.arccos(-2.*val[q]/(1+val[q]**2));
            #
        fie = -np.pi/2.
        #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j)
        #
        ax.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    # real axis
    xi = np.arange(-1, 1, 0.1)
    yi = 0*xi
    plt.plot(xi,yi,color='b',lw=2)
    #
    q = 0
    # text real axis
    while q <= L_val-1:
        xq = (val[q]-1)/(val[q]+1)+0.02;
        ax.text(xq,0.05,0,str(val[q]),color='r',
                 fontweight='normal',fontsize=4)
        q = q+1;
    #
    q = 0
    # text +imag
    while q <= L_val-1:
        #
        xq = (1.j*val[q]-1.)/(1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '+j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='bottom',ha='center',
                         fontweight='normal',fontsize=4)
        elif xp < 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='bottom',ha='right',
                         fontweight='normal',fontsize=4)
        else:
            ax.text(xp,yp,0,sjn,color='b',va='bottom',ha='left',
                         fontweight='normal',fontsize=4)
        q = q+1;
    #
    q = 1
    # text -imag
    while q <= L_val-1:
        #
        xq = (-1.j*val[q]-1.)/(-1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '-j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='top',ha='center',
                         fontweight='normal',fontsize=4)
        elif xp < 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='top',ha='right',
                         fontweight='normal',fontsize=4)
        else:
            ax.text(xp,yp,0,sjn,color='b',va='top',ha='left',
                         fontweight='normal',fontsize=4)
        q = q+1;
    #
#
#####################################################################
#    
def smith3_y(val,ax):
    #
    # Smith chart admittance(y) version for 3d projection
    #
    # val -> values of real and imaginary admittance part circles
    #        example: np.array([0., 0.2, 0.5, 1., 2., 5.])
    #  ax -> axis variable of 3d projection 
    #
    pk = 200
    phi = np.linspace(0,2*np.pi,pk)
    #
    counts = val
    L_val = counts.size
    #
    q = 0
    # circles constant real part
    while q <= L_val-1:
        #
        x0 = -val[q]/(val[q]+1)
        ra = 1./(val[q]+1)
        #
        z = x0 + ra* np.exp(phi*1.j)
        #
        if ra == 1. or ra == 0.5:
            plt.plot(np.real(z),np.imag(z),color='r',
                     lw=2)
        else:
            plt.plot(np.real(z),np.imag(z),color='r',
                     lw=1)
        q = q+1 
    #
    # match circle z -> y g = 1
    z = 0.5 + 0.5*np.exp(phi*1.j);
    plt.plot(np.real(z),np.imag(z),color='m',
             linewidth=1)
    #
    x0 = -1;
    q = 1;
    # circles constant -imaginary part
    while q <= L_val-1:
        #
        y0 = 1./-val[q]
        ra = 1./np.abs(val[q])
        #
        fis = np.pi/2.
        fie = np.pi-np.arccos(-2.*val[q]/(1+val[q]**2))
        #
        if np.abs(val[q]) > 1:
            #
            fie = -np.pi+np.arccos(-2.*val[q]/(1+val[q]**2))
            #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j);
        #
        plt.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    q = 1
    # circles constant +imaginary part
    while q <= L_val-1:
        #
        y0 = 1./val[q]
        ra = 1./np.abs(val[q])
        #
        fis = -np.arccos(2.*val[q]/(1+val[q]**2))
        if np.abs(val[q]) > 1:
            #
            fis = np.pi-np.arccos(-2.*val[q]/(1+val[q]**2));
            #
        fie = -np.pi/2.
        #
        phi = np.linspace(fis,fie,int(pk/np.abs(fie-fis)))
        z = x0+1.j*y0+ra*np.exp(phi*1.j)
        #
        plt.plot(np.real(z),np.imag(z),color='b',
                 lw=1)
        q = q+1
    #
    # real axis
    xi = np.arange(-1, 1.1, 0.1)
    yi = 0*xi
    plt.plot(xi,yi,color='b',lw=2)
    #
    q = 0
    # text real axis
    while q <= L_val-1:
        xq = -(val[q]-1)/(val[q]+1)+0.02;
        ax.text(xq,0.05,0,str(val[q]),color='r',
                 fontweight='normal',fontsize=4)
        q = q+1;
    #
    q = 0
    # text +imag
    while q <= L_val-1:
        #
        xq = -(1.j*val[q]-1.)/(1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '+j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='top',ha='center',
                         fontweight='normal',fontsize=4)
        elif xp < 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='top',ha='right',
                         fontweight='normal',fontsize=4)
        else:
            ax.text(xp,yp,0,sjn,color='b',va='top',ha='left',
                         fontweight='normal',fontsize=4)
        q = q+1;
    #
    q = 1
    # text -imag
    while q <= L_val-1:
        #
        xq = -(-1.j*val[q]-1.)/(-1.j*val[q]+1.)
        xp = np.real(xq)
        yp = np.imag(xq)
        sj = '-j'
        sn = str(val[q])
        sjn = sj+sn
        #
        if xp == 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='bottom',ha='center',
                         fontweight='normal',fontsize=4)
        elif xp < 0:
            #
            ax.text(xp,yp,0,sjn,color='b',va='bottom',ha='right',
                         fontweight='normal',fontsize=4)
        else:
            ax.text(xp,yp,0,sjn,color='b',va='bottom',ha='left',
                         fontweight='normal',fontsize=4)
        q = q+1;
    #
#
#####################################################################










    
    
