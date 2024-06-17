# -*- coding: utf-8 -*-
"""
Spyder-Editor

Dies ist eine temporäre Skriptdatei.
"""
import numpy as np 
#
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
#
[t3sp, t3cs] = t23.t2tot3ce_N(t2sp,t2cs)
print(' ')
print('===== three port parameters common earth =====')
print(' ')
#
print('===== S-parameters =====')
print(' ')
print(t3sp.round(2))
print(' ')
#
print('===== CS-parameters =====')
print(' ')
print(t3cs.round(2))
#
sp_sum_row = np.array([[np.sum(t3sp[0,:])],
                       [np.sum(t3sp[1,:])],
                       [np.sum(t3sp[2,:])]])
#
sp_sum_col = np.array([[np.sum(t3sp[:,0]), np.sum(t3sp[:,1]), np.sum(t3sp[:,2])]])
#
print(' ')
print('===== sum sp-rows =====')
print(' ')
print(sp_sum_row.round(2))
#
print(' ')
print('===== sum sp-columns =====')
print(' ')
print(sp_sum_col.round(2))
#
cs_sum_row = np.array([[np.sum(t3cs[0,:])],
                       [np.sum(t3cs[1,:])],
                       [np.sum(t3cs[2,:])]])
#
cs_sum_col = np.array([[np.sum(t3cs[:,0]), np.sum(t3cs[:,1]), np.sum(t3cs[:,2])]])
#
print(' ')
print('===== sum cs-rows =====')
print(' ')
print(cs_sum_row.round(2))
#
print(' ')
print('===== sum cs-columns =====')
print(' ')
print(cs_sum_col.round(2))
print(' ')
print('========================================================= ')
#
[t3sp, t3cs] = t23.t2tot3ef_N(t2sp,t2cs)
#
print(' ')
print('===== three port parameters earth free =====')
print(' ')
#
print('===== S-parameters =====')
print(' ')
print(t3sp.round(2))
print(' ')
#
print('===== CS-parameters =====')
print(' ')
print(t3cs.round(2))
#
sp_sum_row = np.array([[np.sum(t3sp[0,:])],
                       [np.sum(t3sp[1,:])],
                       [np.sum(t3sp[2,:])]])
#
sp_sum_col = np.array([[np.sum(t3sp[:,0]), np.sum(t3sp[:,1]), np.sum(t3sp[:,2])]])
#
print(' ')
print('===== sum sp-rows =====')
print(' ')
print(sp_sum_row.round(2))
#
print(' ')
print('===== sum sp-columns =====')
print(' ')
print(sp_sum_col.round(2))
#
cs_sum_row = np.array([[np.sum(t3cs[0,:])],
                       [np.sum(t3cs[1,:])],
                       [np.sum(t3cs[2,:])]])
#
cs_sum_col = np.array([[np.sum(t3cs[:,0]), np.sum(t3cs[:,1]), np.sum(t3cs[:,2])]])
#
print(' ')
print('===== sum cs-rows =====')
print(' ')
print(cs_sum_row.round(2))
#
print(' ')
print('===== sum cs-columns =====')
print(' ')
print(cs_sum_col.round(2))
print(' ')
print('========================================================= ')