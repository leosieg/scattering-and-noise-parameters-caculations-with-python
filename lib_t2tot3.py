#
import numpy as np 
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def t2tot3ce_N(t2sp,t2cs):
    #
    # source free two port -> three port, common earth
    #                         signal and noise
    #
    # t2sp -> two port S-parameters emitter/source
    # t2cs -> two port CS-parameters emitter/source
    #
    s2_11 = t2sp[0,0]
    s2_12 = t2sp[0,1]
    s2_21 = t2sp[1,0]
    s2_22 = t2sp[1,1]
    #
    dS = s2_11*s2_22-s2_12*s2_21
    N = 4.-(s2_11+s2_12+s2_21+s2_22)
    #
    s3_11 = (2.*s2_11+1.-dS-s2_12-s2_21)/N
    s3_12 = (2.*s2_12+1.+dS-s2_11-s2_22)/N
    s3_21 = (2.*s2_21+1.+dS-s2_11-s2_22)/N
    s3_22 = (2.*s2_22+1.-dS-s2_12-s2_21)/N
    #
    s3_13 = 1.-(s3_11+s3_12)
    s3_23 = 1.-(s3_21+s3_22)
    s3_31 = 1.-(s3_11+s3_21)
    s3_32 = 1.-(s3_12+s3_22)
    s3_33 = s3_11+s3_12+s3_21+s3_22-1.
    #
    M1 = np.array([[2.-s3_13, -s3_13],
                    [-s3_23, 2.-s3_23],
              [-(1.+s3_33), -(1+s3_33)]])/2.
    #
    M1tc = np.matrix.getH(M1)
    #
    sp3 = np.array([[s3_11, s3_12, s3_13],
                    [s3_21, s3_22, s3_23],
                    [s3_31, s3_32, s3_33]])
    #
    cs3 = np.dot(t2cs,M1tc)
    cs3 = np.dot(M1,cs3)
    #
    return [sp3, cs3]
#
##########################################################################
#
def t2tot3ce_S(t2sp):
    #
    # source free two port -> three port, common earth
    #                         signal only
    #
    # t2sp -> two port S-parameters emitter/source
    #
    s2_11 = t2sp[0,0]
    s2_12 = t2sp[0,1]
    s2_21 = t2sp[1,0]
    s2_22 = t2sp[1,1]
    #
    dS = s2_11*s2_22-s2_12*s2_21
    N = 4.-(s2_11+s2_12+s2_21+s2_22)
    #
    s3_11 = (2.*s2_11+1.-dS-s2_12-s2_21)/N
    s3_12 = (2.*s2_12+1.+dS-s2_11-s2_22)/N
    s3_21 = (2.*s2_21+1.+dS-s2_11-s2_22)/N
    s3_22 = (2.*s2_22+1.-dS-s2_12-s2_21)/N
    #
    s3_13 = 1.-(s3_11+s3_12)
    s3_23 = 1.-(s3_21+s3_22)
    s3_31 = 1.-(s3_11+s3_21)
    s3_32 = 1.-(s3_12+s3_22)
    s3_33 = s3_11+s3_12+s3_21+s3_22-1.
    #
    sp3 = np.array([[s3_11, s3_12, s3_13],
                    [s3_21, s3_22, s3_23],
                    [s3_31, s3_32, s3_33]])
    #
    return sp3
#
##########################################################################
#
def t2tot3ef_N(t2sp,t2cs):
    #
    # source free two port -> three port, earth free
    #                         signal and noise
    #
    # t2sp -> two port S-parameters emitter/source
    # t2cs -> two port CS-parameters emitter/source
    #
    s2_11 = t2sp[0,0]
    s2_12 = t2sp[0,1]
    s2_21 = t2sp[1,0]
    s2_22 = t2sp[1,1]
    #
    dS = s2_11*s2_22-s2_12*s2_21
    N = 4.+s2_11-s2_12-s2_21+s2_22
    #
    s3_11 = (2.*s2_11-1.+dS+s2_12+s2_21)/N
    s3_12 = -(2.*s2_12+1.+dS+s2_11+s2_22)/N
    s3_21 = -(2.*s2_21+1.+dS+s2_11+s2_22)/N
    s3_22 = (2.*s2_22-1.+dS+s2_12+s2_21)/N
    #
    s3_13 = -(1.+s3_11+s3_12)
    s3_23 = -(1.+s3_21+s3_22)
    s3_31 = -(1.+s3_11+s3_21)
    s3_32 = -(1.+s3_12+s3_22)
    s3_33 = 1.+s3_11+s3_12+s3_21+s3_22
    #
    M1 = np.array([[2.+s3_13, -s3_13],
                      [s3_23, -(2.+s3_23)],
                   [s3_33-1., 1.-s3_33]])/2.
    #
    M1tc = np.matrix.getH(M1)
    #
    sp3 = np.array([[s3_11, s3_12, s3_13],
                    [s3_21, s3_22, s3_23],
                    [s3_31, s3_32, s3_33]])
    #
    cs3 = np.dot(t2cs,M1tc)
    cs3 = np.dot(M1,cs3)
    #
    return [sp3, cs3]
#
##########################################################################
#    
def t2tot3ef_S(t2sp):
    #
    # source free two port -> three port, earth free
    #                 signal only
    #
    # t2sp -> two port S-parameters emitter/source
    #
    s2_11 = t2sp[0,0]
    s2_12 = t2sp[0,1]
    s2_21 = t2sp[1,0]
    s2_22 = t2sp[1,1]
    #
    dS = s2_11*s2_22-s2_12*s2_21
    N = 4.+s2_11-s2_12-s2_21+s2_22
    #
    s3_11 = (2.*s2_11-1.+dS+s2_12+s2_21)/N
    s3_12 = -(2.*s2_12+1.+dS+s2_11+s2_22)/N
    s3_21 = -(2.*s2_21+1.+dS+s2_11+s2_22)/N
    s3_22 = (2.*s2_22-1.+dS+s2_12+s2_21)/N
    #
    s3_13 = -(1.+s3_11+s3_12)
    s3_23 = -(1.+s3_21+s3_22)
    s3_31 = -(1.+s3_11+s3_21)
    s3_32 = -(1.+s3_12+s3_22)
    s3_33 = 1.+s3_11+s3_12+s3_21+s3_22
    #
    sp3 = np.array([[s3_11, s3_12, s3_13],
                    [s3_21, s3_22, s3_23],
                    [s3_31, s3_32, s3_33]])
    #
    return sp3
#