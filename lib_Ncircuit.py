#
import numpy as np 
import lib_func.lib_line as ll
import lib_func.lib_conN as cn
import lib_func.lib_spline as spl
import lib_func.lib_t2tot3 as t23
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def lna_circuit_N(freq, dsp):
    #
    # circuit structure low noise amplifier for signal and noise analysis
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
    rcs11 = np.zeros([1, freq_L], dtype = complex)
    rcs12 = np.zeros([1, freq_L], dtype = complex)
    rcs21 = np.zeros([1, freq_L], dtype = complex)
    rcs22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22,cs11,cs12,cs22] = spl.s2spl_N('BFP640F_S.txt',
                                                   'BFP640F_N.txt',freq)
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
        t2cs = np.array([[cs11[q], cs12[q]],
                [np.conj(cs12[q]), cs22[q]]])	   
        #
        # reflection coefficient and noise from emitter inductor
        [gsp,gcs] = ll.gamma_z(jwl*0.5e-9)
        #
        # 3 port transistor matrix with connected inductor, earth free
        [t3sp,t3cs] = t23.t2tot3ce_N(t2sp,t2cs)
        #
        # transistor with emitter inductor
        [rsp,rcs] = cn.Nnpc_xy(t3sp,t3cs,gsp,gcs,3,1,dsp)
        #
        # reflection coefficient and noise from parallel admittance
        # between base and collector
        # here for test purpose only
        # in practice y -> 0
        [gsp,gcs] = ll.gamma_y(1/(1e6/50))
        #
        # 3 port transistor matrix with connected inductor, earth free
        [t3sp,t3cs] = t23.t2tot3ef_N(rsp,rcs)
        #
        # resistor between base and collector
        [rsp,rcs] = cn.Nnpc_xy(t3sp,t3cs,gsp,gcs,3,1,dsp)
        #
        # transistor with Le and parallel admittance
        # change sign, 
        # because different polarity earth free and common earth
        t2sp = np.array([[ rsp[0,0], -rsp[0,1]],
                         [-rsp[1,0],  rsp[1,1]]])
        #
        t2cs = np.array([[ rcs[0,0], -rcs[0,1]],
                         [-rcs[1,0],  rcs[1,1]]])
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
        [rsp,rcs] = cn.Nnpc_xy(lisp,lics,Lisp,Lics,2,1,dsp)
        #
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,t2sp,t2cs,2,1,dsp)
        #
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,Yosp,Yocs,2,1,dsp)
        #
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,Cosp,Cocs,2,1,dsp)
        #
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,losp,locs,2,1,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        rcs11[0,q] = rcs[0,0]
        rcs12[0,q] = rcs[0,1]
        rcs22[0,q] = rcs[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22, rcs11, rcs12, rcs22]
#
#####################################################################
#    
def amp_circuit_N(freq, dsp):
    #
    # circuit structure balanced amplifier for signal and noise analysis
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
    rcs11 = np.zeros([1, freq_L], dtype = complex)
    rcs12 = np.zeros([1, freq_L], dtype = complex)
    rcs21 = np.zeros([1, freq_L], dtype = complex)
    rcs22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22,cs11,cs12,cs22] = spl.s2spl_N('fhc40lg_S.txt','fhc40lg_N.txt',freq)
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
        t2cs = np.array([[cs11[q], cs12[q]],
                [np.conj(cs12[q]), cs22[q]]])	   
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
        [rsp,rcs] = cn.Nnpc_xy(wsp,wcs,t2sp,t2cs,2,1,dsp);
        #
        # Wilkinson-transistor-line
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,lsp,lcs,3,1,dsp)
        #
        #  Wilkinson-line-other-output
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,lsp,lcs,2,1,dsp)
        #
        # Wilkinson-line-other-output-line
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,t2sp,t2cs,3,1,dsp)
        #
        # whole circuit-output Wilkinson
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,wsp,wcs,2,2,dsp)
        #
        # connect inner open ends
        # whole amplifier
        [rsp,rcs] = cn.Nnpc_x(rsp,rcs,2,4,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        rcs11[0,q] = rcs[0,0]
        rcs12[0,q] = rcs[0,1]
        rcs21[0,q] = rcs[1,0]
        rcs22[0,q] = rcs[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22, rcs11, rcs12, rcs22]
#
#####################################################################
#   
def amp_circuit_N02(freq, dsp):
    #
    # circuit structure double balanced amplifier for signal and noise analysis
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
    rcs11 = np.zeros([1, freq_L], dtype = complex)
    rcs12 = np.zeros([1, freq_L], dtype = complex)
    rcs21 = np.zeros([1, freq_L], dtype = complex)
    rcs22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22,cs11,cs12,cs22] = spl.s2spl_N('fhc40lg_S.txt',
                                                   'fhc40lg_N.txt',freq)
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
        t2cs = np.array([[cs11[q], cs12[q]],
                [np.conj(cs12[q]), cs22[q]]])	   
        #
        # Wilkinson splitter 5 GHz
        [wsp,wcs] = ll.wilkinson_s(freq_0,fq)
        #
        # lambda/4-line 50 Ohm
        [lsp,lcs] = ll.line_s(1,l4_0,fq)
        #
        # input circuit
        # lambda/4-line 50 Ohm-transistor
        [ltsp,ltcs] = cn.Nnpc_xy(lsp,lcs,t2sp,t2cs,2,1,dsp)
        #
        # transistor-lambda/4-line 50 Ohm
        [tlsp,tlcs] = cn.Nnpc_xy(t2sp,t2cs,lsp,lcs,2,1,dsp)
        #
        # Wilkinson-Wilkinson
        [rsp,rcs] = cn.Nnpc_xy(wsp,wcs,wsp,wcs,2,1,dsp)
        #
        # 2-Wilkinson-line
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,lsp,lcs,2,1,dsp)
        #
        # 2-Wilkinson-line-Wilkinson
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,wsp,wcs,4,1,dsp);
        #
        # save input circuit
        wlsp = rsp
        wlcs = rcs
        #
        # main amplifier circuit
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,tlsp,tlcs,2,1,dsp)
        #
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,ltsp,ltcs,2,1,dsp)
        #
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,tlsp,tlcs,2,1,dsp)
        #
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,ltsp,ltcs,2,1,dsp)
        #
        # output circuit with saved input circuit
        [rsp,rcs] = cn.Nnpc_xy(rsp,rcs,wlsp,wlcs,2,5,dsp)
        #
        [rsp,rcs] = cn.Nnpc_x(rsp,rcs,2,8,dsp)
        #
        [rsp,rcs] = cn.Nnpc_x(rsp,rcs,2,6,dsp)
        #
        [rsp,rcs] = cn.Nnpc_x(rsp,rcs,2,4,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        rcs11[0,q] = rcs[0,0]
        rcs12[0,q] = rcs[0,1]
        rcs21[0,q] = rcs[1,0]
        rcs22[0,q] = rcs[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22, rcs11, rcs12, rcs22]
#
#####################################################################
#    
def cascode_circuit_N(freq, dsp):
    #
    # circuit structure cascode circuit for signal and noise analysis
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
    rcs11 = np.zeros([1, freq_L], dtype = complex)
    rcs12 = np.zeros([1, freq_L], dtype = complex)
    rcs21 = np.zeros([1, freq_L], dtype = complex)
    rcs22 = np.zeros([1, freq_L], dtype = complex)
    #
    [s11,s12,s21,s22,cs11,cs12,cs22] = spl.s2spl_N('BFP640F_S.txt','BFP640F_N.txt',freq)
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
        t2cs = np.array([[cs11[q], cs12[q]],
                [np.conj(cs12[q]), cs22[q]]])	   
        #
        # 3 port common earth (ce) transistor matrix
        [t3sp,t3cs] = t23.t2tot3ce_N(t2sp,t2cs)
        #
        # short z = 0, reflection coefficient
        [gsp,gcs] = ll.gamma_z(0)
        #
        # transistor short at port 1 == base/gate
        # for common base/gate
        [rsp,rcs] = cn.Nnpc_xy(t3sp,t3cs,gsp,gcs,1,1,dsp)
        #
        # or direct
        #[rsp,rcs] = Nnpc_xy(t3sp,t3cs,-1,0,1,1,dsp);
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
        rcs = np.flip(rcs)
        #
        # cascode connection, collector/drain <---> emitter/source
        [rsp,rcs] = cn.Nnpc_xy(t2sp,t2cs,rsp,rcs,2,1,dsp)
        #
        rsp11[0,q] = rsp[0,0]
        rsp12[0,q] = rsp[0,1]
        rsp21[0,q] = rsp[1,0]
        rsp22[0,q] = rsp[1,1]
        #
        rcs11[0,q] = rcs[0,0]
        rcs12[0,q] = rcs[0,1]
        rcs21[0,q] = rcs[1,0]
        rcs22[0,q] = rcs[1,1]
        #
        q = q + 1
        #
    return [rsp11, rsp12, rsp21, rsp22, rcs11, rcs12, rcs22]        