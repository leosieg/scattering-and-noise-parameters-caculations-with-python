#
import numpy as np 
import matplotlib.pyplot as plt
import lib_func.lib_smith_stab  as sms
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def ncir(noip,val,dsp):
    #
    # plot function noise circles
    #
    # noip -> noise parameters
    #         noip = np.array([[Fmin/dB],
    #                          [Rn/Ohm],
    #                          [CN],
    #                          [Gamma_opt/(re+imj)]])
    #  val -> start point, end point, number of circles-1
    #         example val = np.linspace(0, 3, 7) for 
    #         six circles between Fmin + 3dB
    #         viz. delta = 0.5dB
    #  dsp -> display paramter for Smith chart
    #         'z' -> impedance version, 'y' -> admittance version    
    #
    valq = len(val)
    #
    Fmin = np.real(noip[0,0])
    CN = np.real(noip[2,0]) 
    Gopt = noip[3,0]
    #
    # Schwarz inequality
    Si = 10**(0.1*Fmin)-1.-CN*(1.-(np.abs(Gopt))**2)
    #
    # point density
    pk = 200;
    phi = np.linspace(0,2*np.pi,pk)
    #
    valsmith = np.array([0., 0.2, 0.5, 1., 2., 5.])
    #
    # plot resolution
    plt.rcParams['figure.dpi'] = 600
    #
    # Smith chart z
    if dsp == 'z': 
        sms.smith_z(valsmith)
    #
    # Smith chart y
    if dsp == 'y': 
        sms.smith_y(valsmith)
    #
    plt.axis('off')
    plt.axis('equal') 
    #
    q = 0
    #
    while q <= valq-1:
        #
        n = Fmin*(10**(0.1*val[q])-1.)/CN
        #
        # radius noise circle
        rnc = np.sqrt(n*(1.+n-(np.abs(Gopt))**2))/(1.+n)
        # center noise circle
        cnc = Gopt/(1.+n)
        #
        z = cnc + rnc*np.exp(1j*phi)
        #
        plt.plot(np.real(z), np.imag(z), lw=1, color='k')
        #
        q = q+1
    #
    #Gopt
    xg = np.real(Gopt)
    yg = np.imag(Gopt)
    #
    plt.plot(xg, yg, marker=".", markersize=5)
    #
    plt.text(xg,yg,'$\Gamma_{opt}$',color='k',va='bottom',ha='center',
                             fontweight='bold',fontsize=10) 
    #
    plt.text(-1,1.1,'noise circles, $\Delta F/F_{min}$ = 0.5dB',
             fontsize=8, color='k',fontweight='bold')
    #
    print(' ')
    print('Schwarz inequality Si = ', Si.round(2))
    print(' ')
    #    
#
#####################################################################    
#    
def ncirp(noip,dsp):
    #
    # plot function noise circles and paraboloid
    #
    # noip -> noise parameters
    #         noip = np.array([[Fmin/dB],
    #                          [Rn/Ohm],
    #                          [CN],
    #                          [Gamma_opt/(re+imj)]])
    #  dsp -> display paramter for Smith chart
    #         'z' -> impedance version, 'y' -> admittance version 
    #
    Fmin = np.real(noip[0,0])
    CN = np.real(noip[2,0]) 
    Gopt = noip[3,0]
    #
    valnc = np.linspace(0,3,31)
    Fn = len(valnc)
    #
    valsmith = np.array([0., 0.2, 0.5, 1., 2., 5.])
    #
    val_ncir = np.linspace(0, 3, 7)
    #
    Goptm = np.abs(Gopt)
    Gopta = np.angle(Gopt)
    #
    Fn1 = np.array(np.linspace(0,Fn,Fn+1))
    #
    Fu = Fmin+max(valnc)
    #
    fi = Fmin*(1.+(Fu/Fmin-1.)*np.transpose(np.asmatrix(Fn1))/Fn)
    #
    ni = (np.power(10,0.1*fi)-10**(Fmin*0.1)*np.ones(np.shape(fi)))/CN
    #
    cfi = Goptm/(1.+ni)
    #
    rfi = np.sqrt(np.multiply(ni,(ni+1.-Goptm**2)))/(1+ni)
    #
    #
    Rn = 180.
    #
    fii = np.pi*np.array([np.linspace(-Rn,Rn,int(2*Rn)+1)])/Rn
    #
    fo = np.ones((1,np.size(fii)))*Gopta
    # 
    z1 = np.dot(cfi,np.exp(1.j*fo))
    #
    z2 = np.dot(rfi,np.exp(1.j*fii))
    # 
    Gamma_i = np.add(z1,z2)
    #
    # 
    Gopt_G = Gopt*np.ones(np.shape(Gamma_i))
    # 
    D_G = np.abs(np.add(Gamma_i,-Gopt_G))
    D_G2 = np.multiply(D_G,D_G)
    #
    G_ia = np.abs(Gamma_i)
    G_ia2 = np.multiply(G_ia,G_ia)
    # 
    E = np.ones(np.shape(G_ia))
    # 
    N = np.add(E,-G_ia2)
    # 
    D_GC = CN*np.divide(D_G2,N)
    # 
    Fminl = 10**(0.1*Fmin)*np.ones(np.shape(Gamma_i))
    # 
    F_i = 10*np.log10(np.add(Fminl,D_GC))
    # 
    #################################################
    # plot annotation 
    plt.rcParams.update({"text.usetex": True,
                          "font.family": "Helvetica"})
    #
    #  plot resolution 
    plt.rcParams['figure.dpi'] = 600
    ################################################# 
    # 
    fig = plt.figure()
    # 
    ax = fig.add_subplot(111,projection='3d')
    # 
    ax.plot_surface(np.real(Gamma_i), np.imag(Gamma_i), F_i,cmap=plt.cm.jet)
    # 
    ax.grid(False)
    #
    if dsp == 'z':
        sms.smith3_z(valsmith,ax)
    if dsp == 'y':
        sms.smith3_y(valsmith,ax)    
    # 
    ax.tick_params(axis ='x', labelsize = 6)
    ax.tick_params(axis ='y', labelsize = 6)
    ax.tick_params(axis ='z', labelsize = 6)
    # 
    ax.set_xlabel(r'$x \rightarrow$')
    ax.set_ylabel(r'$y \rightarrow$')
    ax.set_zlabel(r'$F/dB \rightarrow$')
    #
    ax.xaxis.set_ticks(np.arange(-1,1,0.5))
    ax.yaxis.set_ticks(np.arange(-1,1,0.5))
    #
    max_za = int(max(Fmin+val_ncir)+1)
    ax.zaxis.set_ticks(np.arange(0,max_za,0.5))
    #
    pk = 200;
    phi = np.linspace(0,2*np.pi,pk)
    #
    valq = len(val_ncir)
    #  
    q = 0
    # 
    while q <= valq-1:
        #
        n = Fmin*(10**(0.1*val_ncir[q])-1.)/CN
        #
        # radius noise circle
        rnc = np.sqrt(n*(1.+n-(np.abs(Gopt))**2))/(1.+n)
        # center noise circle
        cnc = Gopt/(1.+n)
        #
        cz = cnc + rnc*np.exp(1.j*phi)
        #
        ax.plot(np.real(cz), np.imag(cz), 0, lw=0.75, color='k')
        #
        q = q+1
    #
    #Gopt
    xg = np.real(Gopt)
    yg = np.imag(Gopt)
    #
    ax.plot(xg, yg, 0, marker=".", markersize=5, color = 'k')
    #
    ax.text(xg, yg, 0,'$\Gamma_{opt}$',color ='k',va ='top',ha ='center',
                            fontweight ='bold',fontsize = 10) 
    #
    ax.text(-1.6 , 0, max_za,'noise circles, $\Delta F/F_{min}$ = 0.5dB',
             fontsize = 8, color='k', fontweight='bold')
    #
    return ax
#
#####################################################################       

