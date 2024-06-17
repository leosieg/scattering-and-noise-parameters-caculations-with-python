# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 09:50:13 2024

@author: s_mar
"""
#
import numpy as np 
import matplotlib.pyplot as plt
import sys
#
import lib_func.lib_Ncircuit as cir
import lib_func.lib_smith_stab as sms
import lib_func.lib_noise as noi
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
# single frequency
#freq = np.linspace(1.0, 1.0, 1)*1e9
#dsp = 1
#
# frequency band
freq = np.linspace(0.9, 1.1, 201)*1e9
dsp = 0
#
#
#start = ti.time()
[s11, s12, s21, s22, cs11, cs12, cs22] = cir.lna_circuit_N(freq, dsp)
#ende = ti.time()
#delta_t = ende - start
#print(delta_t)
#
sp = np.array([[s11, s12],[s21, s22]])
cs = np.array([[cs11, cs12],[np.conj(cs12), cs22]])
#
if dsp == 1:
    print('  ')
    print('==========')
    print('=== S ===')
    print(sp.round(2))
    #
    print('  ')
    print('==========')
    print('=== CS ===')
    print(cs.round(2))
#
if dsp == 1:
    sys.exit()
#    
ds = s11*s22-s12*s21
#
muS = (1-np.abs(s22)**2)/(np.abs(s11-np.conj(s22)*ds)+np.abs(s12*s21))
muL = (1-np.abs(s11)**2)/(np.abs(s22-np.conj(s11)*ds)+np.abs(s12*s21)) 
#
# noise parameters
# by means of
# noise wave matrix CS
p = 2.*cs11+2.*(cs22*(1.+np.abs(s11)**2)-2.*np.real(cs12*np.conj(s11)*s21))/np.abs(s21)**2
#
q = 4.*np.abs(cs22*s11-cs12*s21)/np.abs(s21)**2
#
CN = np.real(p+np.sqrt(p**2-q**2))
#
Gopt = 4.*np.conj((cs22*s11-cs12*s21))/(CN*np.abs(s21)**2)
#
Rn = 50.*CN*np.abs(1+Gopt)**2/4.;
#
te_min = 4.*np.real(cs22)/np.abs(s21)**2-CN*np.abs(Gopt)**2
#
Fmin = 10*np.log10(1.+te_min)
#
# noise figure 50 Ohm
F50 = 10*np.log10(10**(Fmin/10.)+CN*np.abs(Gopt)**2);
#
# values at frequency 1e9 Hz
fi = np.where(freq == 1e9)[0][0]
#
sp = np.array([[s11[0,fi], s12[0,fi]],
               [s21[0,fi], s22[0,fi]]])
#
cs = np.array([[cs11[0,fi], cs12[0,fi]],
               [np.conj(cs12[0,fi]), cs22[0,fi]]])
#
noip = np.array([[Fmin[0,fi]],
                [Rn[0,fi]],
                [CN[0,fi]],
                [Gopt[0,fi]]])
#
print('========================')
print('= values at f = 1e9 Hz =')
print(' ')
print('===== S-parameters =====')
print(sp.round(2))
print('========================')
print(' ')
print('===== CS-parameters =====')
print(cs.round(2))
print(' ')
print('=== noise-parameters ===')
print(noip.round(2))
print(' ')
#
s11dB = 20*np.log10(np.abs(s11))
s21dB = 20*np.log10(np.abs(s21))
s22dB = 20*np.log10(np.abs(s22))
#
#######################################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"})
#
# plot resolution 
plt.rcParams['figure.dpi'] = 300
#
plt.figure(1)##########################################
# 
plt.plot(freq/1e9, s21dB[0,:], lw = 4, color = 'b')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{21}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(10, 20)
#plt.savefig('otto.eps', format='eps' ) 
#
plt.figure(2)#########################################
# 
plt.plot(freq/1e9, s11dB[0,:], lw = 4, color = 'r')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{11}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-10, 0)
#
plt.figure(3)#########################################
# 
plt.plot(freq/1e9, s22dB[0,:], lw = 4, color = 'r')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{22}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-40, 0)
#
plt.figure(4)######################################### 
#
plt.plot(freq/1e9, muS[0,:], lw = 4, color ='r')
plt.plot(freq/1e9, muL[0,:], lw = 4, color ='b')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$\mu \rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(0, 5)
#
plt.legend(['$\mu_S$', '$\mu_L$'],loc = 'lower right', frameon = False)
#
plt.figure(5)#########################################
#
sms.ssks(sp,'z')
#
plt.figure(6)#########################################
#
sms.sskl(sp,'y')
#
plt.figure(7)#########################################
#
plt.plot(freq/1e9, Fmin[0,:], lw = 4, color = 'r')
plt.plot(freq/1e9, F50[0,:], lw = 1.5, color = 'b')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'F/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(0.6, 0.7)
#
plt.legend(['$F_{min}$', '$F_{\Gamma_S=0}$'],loc = 'lower right', frameon = False)
#   
plt.figure(8)#########################################
#
val = np.linspace(0, 3, 7)
noi.ncir(noip, val, 'z')
#
plt.figure(9)#########################################
#
noi.ncirp(noip, 'y')
#


