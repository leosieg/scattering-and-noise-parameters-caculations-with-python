#
import numpy as np
from scipy.interpolate import CubicSpline 
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def s2spl_S(data_S,freq):
    #
    # spline approximation for data in tables *.s2p,
    #              signal only
    #
    # data_S -> '*.txt' file with two port data for frequency and 
    #           S-parameters in mag/deg !!! format
    #           example: 'BF640F_S.txt'
    #   freq -> frequency/Hz
    #          
    dat = data_S
    #
    file = './transistors/' + dat
    #
    Sp = np.loadtxt(file)
    #
    sp = np.zeros([len(freq),8])
    #
    freq_Sp = Sp[:,0]*1.e9
    #
    q = 1
    #
    while q <= 8:
         fc = CubicSpline(freq_Sp, Sp[:, q], bc_type='natural')
         sp[:, q-1] = fc(freq)
         q = q+1
    #
    s11 = sp[:,0]*np.exp(-sp[:,1]*np.pi/180.j)
    s21 = sp[:,2]*np.exp(-sp[:,3]*np.pi/180.j)
    s12 = sp[:,4]*np.exp(-sp[:,5]*np.pi/180.j)
    s22 = sp[:,6]*np.exp(-sp[:,7]*np.pi/180.j)
    #
    return [s11, s12, s21, s22] 
#
##########################################################################
#
def s2spl_N(data_S,data_N,freq):
    #
    # spline approximation for data in tables *.s2p,
    #              signal and noise only
    #
    # data_S -> '*.txt' file with 50 Ohm two port data for frequency and 
    #           S-parameters in mag/deg !!! format
    #           example: 'BF640F_S.txt'
    # data_N -> '*.txt' file with 50 Ohm two port data for frequency and 
    #           noise parameters in dB, mag/deg, normalized noise resistance !!! format
    #           example: 'BF640F_N.txt'
    #   freq -> frequency/Hz
    #    
    dat = data_S
    #
    file = './transistors/' + dat
    #
    Sp = np.loadtxt(file)
    #
    sp = np.zeros([len(freq),8])
    #
    freq_Sp = Sp[:,0]*1.e9
    #
    q = 1
    #
    while q <= 8:
         fc = CubicSpline(freq_Sp, Sp[:, q], bc_type='natural')
         sp[:, q-1] = fc(freq)
         q = q+1
    #
    s11 = sp[:,0]*np.exp(-sp[:,1]*np.pi/180.j)
    s21 = sp[:,2]*np.exp(-sp[:,3]*np.pi/180.j)
    s12 = sp[:,4]*np.exp(-sp[:,5]*np.pi/180.j)
    s22 = sp[:,6]*np.exp(-sp[:,7]*np.pi/180.j)
    #
    dat = data_N
    #
    file = './transistors/' + dat
    #
    Noip = np.loadtxt(file)  
    #
    noip = np.zeros([len(freq),4])
    #
    freq_Np = Noip[:,0]*1.e9
    #
    q = 1
    #
    while q <= 4:
         fc = CubicSpline(freq_Np, Noip[:, q], bc_type='natural')
         noip[:, q-1] = fc(freq)
         q = q+1
    #
    Fmin = noip[:,0]
    Gopt = noip[:,1]*np.exp(-noip[:,2]*np.pi/180.j)
    Rn = 50.*noip[:,3]
    CN = 4.*Rn/(50.*np.abs(1+Gopt)**2)
    #
    te_min = 10.**(0.1*Fmin)-1.
    #
    cs11 = (te_min*(np.abs(s11)**2-1.)+CN*np.abs(1.-s11*Gopt)**2)/4.;
    cs12 = np.conj(s21)*(s11*te_min-CN*np.conj(Gopt)*(1.-s11*Gopt))/4.;
    cs22 = (np.abs(s21)**2*(te_min+CN*np.abs(Gopt)**2))/4.;
    #
    return [s11,s12,s21,s22,cs11,cs12,cs22]