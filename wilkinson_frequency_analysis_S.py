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
import lib_func.lib_Scircuit as cir
import lib_func.lib_smith_stab as sms
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
# single frequency
#freq = np.linspace(1.0, 1.0, 1)*1e9
#freq = [1e9]
#dsp = 1
#
# frequency band
freq = np.linspace(0, 2, 501)*1e9
dsp = 0
#
#
#start = ti.time()
[s11, s12, s13, s22, s23, s33] = cir.wilkinson_s(1e9, freq, dsp)
#ende = ti.time()
#delta_t = ende - start
#print(delta_t)
#
#
sp = np.array([[s11, s12, s13],
               [s12, s22, s23],
               [s13, s23, s33]])
#
if dsp == 1:
    print('  ')
    print('==========')
    print('=== S ===')
    print(sp.round(2))
    print('  ')
#
if dsp == 1:
    sys.exit()
#
# values at frequency 1e9 Hz
fi = np.where(freq == 1e9)[0][0]
#
sp = np.array([[s11[0,fi], s12[0,fi], s13[0,fi]],
               [s12[0,fi], s22[0,fi], s23[0,fi]],
               [s13[0,fi], s23[0,fi], s33[0,fi]]])
#
print('========================')
print('= values at f = 1e9 Hz =')
print(' ')
print('===== S-parameters =====')
print(sp.round(2))
print('========================')
print(' ')
#
#
s11dB = 20*np.log10(np.abs(s11))
s12dB = 20*np.log10(np.abs(s12))
s13dB = 20*np.log10(np.abs(s13))
#
s22dB = 20*np.log10(np.abs(s22))
s23dB = 20*np.log10(np.abs(s23))
#
s33dB = 20*np.log10(np.abs(s33))
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
plt.plot(freq/1e9, s11dB[0,:], lw = 4, color = 'b')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{11}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-50, 0)
#plt.savefig('otto.eps', format='eps' ) 
#
plt.figure(2)#########################################
# 
plt.plot(freq/1e9, s12dB[0,:], lw = 4, color = 'r')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{12}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-5, 0)
#
plt.figure(3)#########################################
# 
plt.plot(freq/1e9, s13dB[0,:], lw = 4, color = 'r')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{13}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-5, 0)
#
plt.figure(4)#########################################
# 
plt.plot(freq/1e9, s22dB[0,:], lw = 4, color = 'r')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{22}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-50, 0)
#
plt.figure(5)#########################################
# 
plt.plot(freq/1e9, s23dB[0,:], lw = 4, color = 'r')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{23}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-50, 0)
#
plt.figure(6)#########################################
# 
plt.plot(freq/1e9, s33dB[0,:], lw = 4, color = 'r')
#
plt.xlabel(r'frequency/GHz $\rightarrow$', fontsize = 12)
plt.ylabel(r'$|s_{33}|^2$/dB $\rightarrow$', fontsize = 12)
plt.grid()
plt.xlim(min(freq)/1e9, max(freq)/1e9)
plt.ylim(-50, 0)
#

