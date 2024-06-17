#
import numpy as np 
import lib_func.lib_line  as ll
import lib_func.lib_conS  as cs
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def line_s(zl,ll,freq):
    #
    # mismatched homogeneous series line
    #
    #   zl -> normalized chracteristic impedance/50 Ohm
    #   ll -> line length/m
    # freq -> frequency/Hz
    #
    c0 = 2.99792458e8
    bl = 2*np.pi*freq*ll/c0
    g = (zl-1.)/(zl+1.)
    #
    n = 1-g*g*np.exp(-2.j*bl)
    #
    rho = g*(1-np.exp(-2.j*bl))/n
    tau = np.exp(-1.j*bl)*(1.-g*g)/n
    #
    sp = np.array([[rho, tau],
                   [tau, rho]])
    #
    sptc = np.matrix.getH(sp)
    sptc = np.dot(sp,sptc)
    #
    cs = 0.*np.add(np.eye(2, dtype=float),-sptc)/4.
    #    
    return [sp, cs] 
#
##########################################################################
#
def line_p(zl,ll,freq,en):
    #
    #  parallel stub line
    #
    #   zl -> normalized chracteristic impedance/50 Ohm
    #   ll -> line length/m
    # freq -> frequency/Hz
    #   en -> end indicator, 'o' -> open, 's' -> short
    #
    c0 = 2.99792458e8
    #
    if en == 'o':
        y = 1.j*np.tan(2.*np.pi*freq*ll/c0)/zl
    #
    if en == 's':
        y = -1.j/(zl*np.tan(2.*np.pi*freq*ll/c0))
    #
    sp = np.array([[-y, 2.],
                   [2., -y]])/(2.+y)
    #
    sptc = np.matrix.getH(sp)
    sptc = np.dot(sp,sptc)
    #
    cs = 0.*np.add(np.eye(2, dtype=float),-sptc)/4.
    #    
    return [sp, cs]
#
##########################################################################
#    
def gamma_z(z):
    #
    # reflection coefficient for z
    #
    # z -> normalized impedance
    #
    gsp = np.zeros((1,1), dtype = complex)
    #
    gsp[0,0] = (z-1.)/(z+1.)
    #
    gsptc = np.matrix.getH(gsp)
    gsptc = np.dot(gsp,gsptc)
    #
    gcs = np.add(np.eye(1, dtype=float),-gsptc)/4.
    #
    return [gsp, gcs]	
#
##########################################################################
#     
def gamma_y(y):
    #
    # reflection coefficient for y
    #
    # y -> normalized admittance
    #
    gsp = np.zeros((1,1), dtype = complex)
    #
    gsp[0,0] = -(y-1.)/(y+1.)
    #
    gsptc = np.matrix.getH(gsp)
    gsptc = np.dot(gsp,gsptc)
    #
    gcs = np.add(np.eye(1, dtype=float),-gsptc)/4.
    #
    return [gsp, gcs]	
#
##########################################################################
#   
def y_parallel(y):
    #
    # parallel admittance
    #
    # y -> normalized admittance
    #
    ysp = np.array([[-y, 2.],
                    [2., -y]])/(2.+y)
    # 
    sptc = np.matrix.getH(ysp)
    sptc = np.dot(ysp,sptc)
    #
    ycs = np.add(np.eye(2,dtype=float),-sptc)/4.
    #
    return [ysp, ycs]
#
##########################################################################
#  	
def z_series(z):
    #
    # series impedance
    #
    # z -> normalized impedance
    #
    zsp = np.array([[z, 2.],
                    [2., z]])/(2.+z)
    # 
    sptc = np.matrix.getH(zsp)
    sptc = np.dot(zsp,sptc)
    #
    zcs = np.add(np.eye(2,dtype=float),-sptc)/4.
    #
    return [zsp, zcs]
#
##########################################################################
#  	
def jomegalc(freq):
    #
    # factor for normalizing inductance or capacitance
    #
    # freq -> frequency/Hz
    #    
    jwl = 2.j*np.pi*freq/50.
    jwc = 2.j*np.pi*freq*50.
    #
    return [jwl, jwc]
#
##########################################################################
#     
def p3sp():
    # 
    # three port splitter, common earth
    #
    sp = np.array([[-1., 2., 2.],
                   [2., -1., 2.],
	               [2., 2., -1.]])/3.
    #
    cs = 0.*sp	 
    #
    return [sp, cs]
#
##########################################################################
# 
def wilkinson_s(freq_0, freq):
    #
    # 50 Ohm Wilkinson splitter
    #
    # freq_0 -> design frequency/Hz
    #   freq -> current frequency/Hz
    #
    c0 = 2.99792458e8;
    # lambda/4 for freq_0
    l4_0 = c0/(4.*freq_0)
    #
    [spx,csx] = ll.p3sp()
    psp = spx
    #
    # lambda/4-line ZL/50 = sqrt(2)
    # at the design frequency
    [spl, csl] = ll.line_s(np.sqrt(2.),l4_0,freq)
    #
    rsp = cs.Snpc_xy(spx,spl,2,1,0)
    #
    rsp = cs.Snpc_xy(rsp,spl,2,1,0)
    #
    rsp = cs.Snpc_xy(rsp,psp,2,1,0)
    #
    rsp = cs.Snpc_xy(rsp,psp,2,1,0)
    #
    # Wilkinson resistor 100 ohm
    [spz,csz] = ll.z_series(2.)
    #
    rsp = cs.Snpc_xy(rsp,spz,3,1,0)
    #
    # S-paramters whole Wilkinson splitter
    rsp = cs.Snpc_x(rsp,3,5,0)
    #
    rsptc = np.matrix.getH(rsp)
    #
    rcs = np.dot(rsp,rsptc)
    #
    # noise wave matrix whole Wilkinson splitter (passive n-port)
    rcs = np.add(np.eye(3, dtype=float),-rcs)/4.
    #
    return [rsp, rcs]