# -*- coding: utf-8 -*-
"""
Spyder-Editor

Dies ist eine temporäre Skriptdatei.
"""
import numpy as np 
#
import lib_func.lib_line  as ll
import lib_func.lib_conN  as cn
import lib_func.lib_spline  as spl
import lib_func.lib_t2tot3  as t23
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
# frequency
freq = np.linspace(1.0, 1.0, 1)*1e9
#
# spline with data for transistor common emitter
[s11,s12,s21,s22,cs11,cs12,cs22] = spl.s2spl_N('BFP640F_S.txt','BFP640F_N.txt',freq)
#
print(' ')
print('===== two port parameters =====')
print(' ')
# 2 port transistor matrices common emitter
t2sp = np.array([[s11[0], s12[0]],
                 [s21[0], s22[0]]])
#
t2cs = np.array([[cs11[0], cs12[0]],
        [np.conj(cs12[0]), cs22[0]]])
#
print('===== S-parameters =====')
print(' ')
print(t2sp.round(2))
print(' ')
#
print('===== CS-parameters =====')
print(' ')
print(t2cs.round(2))
print(' ')
#
# 3 port transistor matrix, common earth
[t3sp, t3cs] = t23.t2tot3ce_N(t2sp,t2cs)
#
# normalized reactance/suszeptance
[jwl, jwc] = ll.jomegalc(freq)
#
# reflection coefficient and noise from emitter inductor
[gsp, gcs] = ll.gamma_z(jwl*0.5e-9)
#
# transistor with emitter inductor
[rsp, rcs] = cn.Nnpc_xy(t3sp,t3cs,gsp,gcs,3,1,0)
#
print(' ')
print('===== inductor at emitter =====')
print(' ')
print('===== S-parameters =====')
print(' ')
print(rsp.round(2))
print(' ')
#
print('===== CS-parameters =====')
print(' ')
print(rcs.round(2))
print(' ')
#
print(' ')
print('===== transistor with emitter inductor as earth free 3 port !!! =====')
print(' ')
#3 port transistor matrix with connected inductor, earth free
[t3sp, t3cs] = t23.t2tot3ef_N(rsp,rcs)
#
print(' ')
print('===== S-parameters =====')
print(' ')
print(t3sp.round(2))
print(' ')
#
print('===== CS-parameters =====')
print(' ')
print(t3cs.round(2))
print(' ')
#
# reflection coefficient and noise from parallel admittance
# between base and collector
#[gsp,gcs] = ll.gamma_y(1/(300/50))
[gsp, gcs] = ll.gamma_y(1/(300/50))
#
print(' ')
print('===== resistor between base/gate and collector/drain =====')
print('===== as an 1 port with reflection coefficient Gamma =====')
print(' ')
#
[rsp, rcs] = cn.Nnpc_xy(t3sp,t3cs,gsp,gcs,3,1,1)
#
print(' ')
print('===== S-parameters =====')
print(' ')
print(rsp.round(2))
print(' ')
#
print('===== CS-parameters =====')
print(' ')
print(rcs.round(2))
print(' ')
#
# change sign transmission,
# because different polarity earth free and common earth
print(' ')
print('===================== change sign ========================')
print(' ')
#
rsp = np.array([[rsp[0,0], -rsp[0,1]],
              [-rsp[1,0], rsp[1,1]]])
#
rcs = np.array([[rcs[0,0], -rcs[0,1]],
                [-rcs[1,0], rcs[1,1]]]) 
#
print(' ')
print('===== S-parameters =====')
print(' ')
print(rsp.round(2))
print(' ')
#
print('===== CS-parameters =====')
print(' ')
print(rcs.round(2))
print(' ')
#
s11 = rsp[0,0]
s21 = rsp[1,0]
#
cs11 = rcs[0,0]
cs12 = rcs[0,1]
cs22 = rcs[1,1]
#
# noise parameters
# by means of
# noise wave matrix CS
p = 2.*cs11+2.*(cs22*(1+np.abs(s11)**2)-2.*np.real(cs12*np.conj(s11)*s21))/np.abs(s21)**2
#
q = 4.*np.abs(cs22*s11-cs12*s21)/np.abs(s21)**2
#
CN = np.real(p+np.sqrt(p**2-q**2))
#
Gopt = 4.*np.conj((cs22*s11-cs12*s21))/(CN*np.abs(s21)**2)
#
Rn = 50.*CN*np.abs(1+Gopt)**2/4.
#
te_min = 4.*np.real(cs22)/np.abs(s21)**2-CN*np.abs(Gopt)**2
#
Fmin = 10*np.log10(1+te_min)
#
noip = np.array([[Fmin],
                 [Rn],
                 [CN],
                 [Gopt]])
#
print(' ')
print('===== noise-parameters =====')
print(' ')
print(noip.round(2))
print(' ')
print('============================================================')
#
