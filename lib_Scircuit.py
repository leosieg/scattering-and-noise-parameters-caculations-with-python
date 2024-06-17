#
import numpy as np 
import lib_func.lib_line as ll
import lib_func.lib_conS as cs
import lib_func.lib_spline as spl
import lib_func.lib_t2tot3 as t23
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def lna_circuit_S(freq, dsp):
    #
    # circuit structure low noise amplifier for signal analysis only
    #
    # freq -> frequency/Hz
    # dsp -> display parameter to show the port connection numbering,
    #        should be used for single frequency analysis only
    #
    freq_L = len(freq)
    #
    # matrix adjustments
    rsp11 = np.zeros([1, freq_L], dtype = complex)
    rsp12 = np.zeros([1, freq_L], dtype = complex)
    rsp21 = np.zeros([1, freq_L], dtype = complex)
    rsp22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22] = spl.s2spl_S('BFP640F_S.txt',freq)
    #
    c0 = 2.99792458e8;
    # lambda/4 for freq = 1e9
    l4_0 = c0/(4*1e9);
    #
    q = 0
    #
    while q <= freq_L-1:
        #
        fq = freq[q]
        #
        # nomalized inductance, capacitance 
        [jwl,jwc] = ll.jomegalc(fq)    
        #
        # transistor
        t2sp = np.array([[s11[q], s12[q]],
                          [s21[q], s22[q]]])   
        #
        # reflection coefficient and noise from emitter inductor
        [gsp,gcs] = ll.gamma_z(jwl*0.5e-9)
        #
        # 3 port transistor matrix with connected inductor, earth free
        t3sp = t23.t2tot3ce_S(t2sp)
        #
        # transistor with emitter inductor
        rsp = cs.Snpc_xy(t3sp,gsp,3,1,dsp)
        #
        # reflection coefficient and noise from parallel admittance
        # between base and collector
        # for test purpose only
        # in practice y -> 0
        [gsp,gcs] = ll.gamma_y(1/(1e6/50))
        #
        # 3 port transistor matrix with connected inductor, earth free
        t3sp = t23.t2tot3ef_S(rsp)
        #
        # resistor between base and collector
        rsp = cs.Snpc_xy(t3sp,gsp,3,1,dsp)
        #
        # transistor with Le and parallel admittance
        # change sign, 
        # because different polarity earth free and common earth
        t2sp = np.array([[ rsp[0,0], -rsp[0,1]],
                         [-rsp[1,0],  rsp[1,1]]])
        #
        # input lambda/4-line 78 Ohm
        [lisp,lics] = ll.line_s(78/50,l4_0,fq)
        #
        # input parallel L 
        [Lisp,Lics] = ll.y_parallel(1/(jwl*55.5e-9))
        #
        # output parallel y
        [Yosp,Yocs] = ll.y_parallel(1/(30/50+jwl*2e-9))
        #
        # output parallel C
        [Cosp,Cocs] = ll.y_parallel(jwc*0.5e-12)
        #
        # ouput lambda/4-line 38 Ohm
        [losp,locs] = ll.line_s(38/50,l4_0,fq)
        #
        # amplifier circuit
        #
        rsp = cs.Snpc_xy(lisp,Lisp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,t2sp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,Yosp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,Cosp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,losp,2,1,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22]
#
#####################################################################
#    
def amp_circuit_S(freq, dsp):
    #
    # circuit structure balanced amplifier for signal analysis only
    #
    # freq -> frequency/Hz
    # dsp -> display parameter to show the port connection numbering,
    #        should be used for single frequency analysis only
    #
    freq_L = len(freq)
    #
    # matrix adjustments
    rsp11 = np.zeros([1, freq_L], dtype = complex)
    rsp12 = np.zeros([1, freq_L], dtype = complex)
    rsp21 = np.zeros([1, freq_L], dtype = complex)
    rsp22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22] = spl.s2spl_S('fhc40lg_S.txt',freq)
    #
    freq_0 = 5e9;
    #
    c0 = 2.99792458e8;
    # lambda/4 for freq_0
    l4_0 = c0/(4*freq_0);
    #
    q = 0
    #
    while q <= freq_L-1:
        #
        fq = freq[q]
        #
        # transistor
        t2sp = np.array([[s11[q], s12[q]],
                         [s21[q], s22[q]]])
        #
        # Wilkinson splitter 5 GHz
        [wsp,wcs] = ll.wilkinson_s(freq_0,fq)
        #
        # lambda/4-line 50 Ohm
        [lsp,lcs] = ll.line_s(1,l4_0,fq)
        #
        # amplifier circuit
        #
        # input Wilkinson-transistor
        rsp = cs.Snpc_xy(wsp,t2sp,2,1,dsp);
        #
        # Wilkinson-transistor-line
        rsp = cs.Snpc_xy(rsp,lsp,3,1,dsp)
        #
        #  Wilkinson-line-other-output
        rsp = cs.Snpc_xy(rsp,lsp,2,1,dsp)
        #
        # Wilkinson-line-other-output-line
        rsp = cs.Snpc_xy(rsp,t2sp,3,1,dsp)
        #
        # whole circuit-output Wilkinson
        rsp = cs.Snpc_xy(rsp,wsp,2,2,dsp)
        #
        # connect inner open ends
        # whole amplifier
        rsp = cs.Snpc_x(rsp,2,4,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22]
#
#####################################################################
#   
def amp_circuit_S02(freq, dsp):
    #
    # circuit structure double balanced amplifier for signal analysis only
    #
    # freq -> frequency/Hz
    # dsp -> display parameter to show the port connection numbering,
    #        should be used for single frequency analysis only
    #
    freq_L = len(freq)
    #
    # matrix adjustments
    rsp11 = np.zeros([1, freq_L], dtype = complex)
    rsp12 = np.zeros([1, freq_L], dtype = complex)
    rsp21 = np.zeros([1, freq_L], dtype = complex)
    rsp22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22] = spl.s2spl_S('fhc40lg_S.txt',freq)
    #
    # freq = freq_0
    freq_0 = 5e9;
    #
    c0 = 2.99792458e8;
    # lambda/4 for freq_0
    l4_0 = c0/(4*freq_0);
    #
    q = 0
    #
    while q <= freq_L-1:
        #
        fq = freq[q]
        #
        # transistor
        t2sp = np.array([[s11[q], s12[q]],
                         [s21[q], s22[q]]])
        #
        # Wilkinson splitter 5 GHz
        [wsp,wcs] = ll.wilkinson_s(freq_0,fq)
        #
        # lambda/4-line 50 Ohm
        [lsp,lcs] = ll.line_s(1,l4_0,fq)
        #
        # input circuit
        # lambda/4-line 50 Ohm-transistor
        ltsp = cs.Snpc_xy(lsp,t2sp,2,1,dsp)
        #
        # transistor-lambda/4-line 50 Ohm
        tlsp = cs.Snpc_xy(t2sp,lsp,2,1,dsp)
        #
        # Wilkinson-Wilkinson
        rsp = cs.Snpc_xy(wsp,wsp,2,1,dsp)
        #
        # 2-Wilkinson-line
        rsp= cs.Snpc_xy(rsp,lsp,2,1,dsp)
        #
        # 2-Wilkinson-line-Wilkinson
        rsp = cs.Snpc_xy(rsp,wsp,4,1,dsp);
        #
        # save input circuit
        wlsp = rsp
        #
        # main amplifier circuit
        rsp = cs.Snpc_xy(rsp,tlsp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,ltsp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,tlsp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,ltsp,2,1,dsp)
        #
        # output circuit with saved input circuit
        rsp = cs.Snpc_xy(rsp,wlsp,2,5,dsp)
        #
        rsp = cs.Snpc_x(rsp,2,8,dsp)
        #
        rsp = cs.Snpc_x(rsp,2,6,dsp)
        #
        rsp = cs.Snpc_x(rsp,2,4,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22]
#
#####################################################################
#    
def cascode_circuit_S(freq, dsp):
    #
    # circuit structure cascode circuit for signal analysis only 
    #
    # freq -> frequency/Hz
    # dsp -> display parameter to show the port connection numbering,
    #        should be used for single frequency analysis only
    #
    freq_L = len(freq)
    #
    # matrix adjustments
    rsp11 = np.zeros([1, freq_L], dtype = complex)
    rsp12 = np.zeros([1, freq_L], dtype = complex)
    rsp21 = np.zeros([1, freq_L], dtype = complex)
    rsp22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22] = spl.s2spl_S('BFP640F_S.txt',freq)
    #
    q = 0
    #
    while q <= freq_L-1:
        #
        fq = freq[q]
        #
        # transistor
        t2sp = np.array([[s11[q], s12[q]],
                         [s21[q], s22[q]]])
        #
        # 3 port common earth (ce) transistor matrix
        t3sp = t23.t2tot3ce_S(t2sp)
        #
        # short z = 0, reflection coefficient
        [gsp,gcs] = ll.gamma_z(0)
        #
        # transistor short at port 1 == base/gate
        # for common base/gate
        rsp = cs.Snpc_xy(t3sp,gsp,1,1,dsp)
        #
        # or direct
        # rsp = cs.Snpc_xy(t3sp,-1,1,1,dsp);
        #
        # because in result after connection 
        # old port 2 collector = new port 1 
        # and 
        # old port 3 emitter = new port 2
        # but it must be
        # port 1 emitter, port 2 collector
        # for common base
        # we must use flip to get a matrix from
        # [a b    to      [d c
        #  c d]            b a]
        #
        # common base/gate
        rsp = np.flip(rsp)
        #
        # cascode connection, collector/drain <---> emitter/source
        rsp = cs.Snpc_xy(t2sp,rsp,2,1,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22]
#
#####################################################################
# 
def filter_circuit_S(freq, dsp):
    #
    # circuit structure filter circuit for signal analysis only
    #
    # freq -> frequency/Hz
    # dsp -> display parameter to show the port connection numbering,
    #        should be used for single frequency analysis only
    #
    freq_L = len(freq)
    #
    # matrix adjustments
    rsp11 = np.zeros([1, freq_L], dtype = complex)
    rsp12 = np.zeros([1, freq_L], dtype = complex)
    rsp21 = np.zeros([1, freq_L], dtype = complex)
    rsp22 = np.zeros([1, freq_L], dtype = complex)
    #
    q = 0
    #
    while q <= freq_L-1:
        #
        fq = freq[q]
        #
        # nomalized inductance, capacitance 
        [jwl,jwc] = ll.jomegalc(fq)    
        #
        # each resonant circuit
        [sp1,cs1] = ll.y_parallel(jwc*10.659e-12+1/(jwl*2.4755e-9))
        #
        [sp2,cs2] = ll.z_series(jwl*26.5991e-9+1/(jwc*991.9785e-15))
        #
        [sp3,cs3] = ll.y_parallel(jwc*17.2369e-12+1/(jwl*1.5308e-9))
        #
        [sp4,cs4] = ll.z_series(jwl*26.5991e-9+1/(jwc*991.9785e-15))
        #
        [sp5,cs5] = ll.y_parallel(jwc*10.659e-12+1/(jwl*2.4755e-9))
        #
        rsp = cs.Snpc_xy(sp1,sp2,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,sp3,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,sp4,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,sp5,2,1,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22]    
#
#####################################################################
# 
def wilkinson_s(freq_0, freq, dsp):
    #
    freq_L = len(freq)
    #
    # matrix adjustments
    rsp11 = np.zeros([1, freq_L], dtype = complex)
    rsp12 = np.zeros([1, freq_L], dtype = complex)
    rsp13 = np.zeros([1, freq_L], dtype = complex)
    rsp22 = np.zeros([1, freq_L], dtype = complex)
    rsp23 = np.zeros([1, freq_L], dtype = complex)
    rsp33 = np.zeros([1, freq_L], dtype = complex)
    #
    c0 = 2.99792458e8;
    # lambda/4 for freq_0
    l4_0 = c0/(4.*freq_0)
    #
    [spx,csx] = ll.p3sp()
    psp = spx
    #
    q = 0
    #
    while q <= freq_L-1:
        #
        fq = freq[q]
        #
        # lambda/4-line ZL/50 = sqrt(2)
        # at the design frequency
        [spl, csl] = ll.line_s(np.sqrt(2),l4_0,fq)
        #
        rsp = cs.Snpc_xy(spx,spl,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,spl,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,psp,2,1,dsp)
        #
        rsp = cs.Snpc_xy(rsp,psp,2,1,dsp)
        #
        # Wilkinson resistor 100 ohm
        [spz,csz] = ll.z_series(2.)
        #
        rsp = cs.Snpc_xy(rsp,spz,3,1,dsp)
        #
        # S-paramters whole Wilkinson splitter
        rsp = cs.Snpc_x(rsp,3,5,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp13[0,q] = rsp[0,2]
        rsp22[0,q] = rsp[1,1]
        rsp23[0,q] = rsp[1,2]
        rsp33[0,q] = rsp[2,2]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp13, rsp22, rsp23, rsp33]     