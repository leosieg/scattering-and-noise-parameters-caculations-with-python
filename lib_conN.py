#
import numpy as np
import sys
#
# Siegfried Martius
# siegfried.martius@fau.de
# www.lhft.de
#
def Nnpc_xy(Sx,CSx,Sy,CSy,P1x,P2y,dsp):
    #
    # port connection x-port -> y-port, signal and noise
    #
    #  Sx -> S-matrix x-port
    # CSx -> CS-matrix x-port
    #  Sy -> S-matrix y-port
    # CSy -> CS-matrix y-port
    # P1x -> port number from x-port for connection with P2y
    # P2y -> port number from y-port for connection with P1x
    # dsp -> display paramter to show the port numbers of connection
    #        before and after
    #
    Px = len(Sx)
    Py = len(Sy)
    #
    S_1 = np.zeros((Px-1+Py-1,Px-1+Py-1), dtype = complex)
    CS_1 = np.zeros((Px-1+Py-1,Px-1+Py-1), dtype = complex)
    #
    S_2 = np.zeros((Px-1+Py-1,2), dtype = complex);
    CS_2 = np.zeros((Px-1+Py-1,2), dtype = complex);
    #
    S_3 = np.zeros((2,Px-1+Py-1), dtype = complex);
    #
    Nx = np.array(list(range(1, Px+1)))
    Nx0 = Nx
    #
    Nx = np.delete(Nx, P1x-1)
    #
    Ny = np.array(list(range(1, Py+1)))
    Ny0 = Ny
    #
    Ny = np.delete(Ny, P2y-1)
    #
    Nuv = np.append(Nx, Ny)
    #
    Pc = [P1x, P2y]
    #
    Nuv = np.append(Nuv, Pc)
    #
    if dsp == 1:
       N_Tor = len(Nuv)-2
       N_Tor_n = np.array(list(range(1,N_Tor+1)))
       #
       print('==========================================')
       print('===== connection x-port <--> y-port ======')
       print('port number x: ', P1x, '; port number y: ', P2y, ';')
       print('port numbers start: ', Nx0, ' |*| ', Ny0)
       print('port numbers  open: ', Nx, ' |*| ', Ny)
       print('port numbers   new: ', N_Tor_n)
       print('==========================================')   
    #
    m = 0
    while m <= Px-2:
        n = 0
        while n <= Px-2:
            S_1[m,n] = Sx[Nuv[m]-1,Nuv[n]-1]
            CS_1[m,n] = CSx[Nuv[m]-1,Nuv[n]-1]
            n = n+1
        m = m+1
    #
    m = Px-1
    while m <= Px-2+Py-1:
        n = Px-1
        while n <= Px-2+Py-1:
            S_1[m,n] = Sy[Nuv[m]-1,Nuv[n]-1]
            CS_1[m,n] = CSy[Nuv[m]-1,Nuv[n]-1]
            n = n+1
        m = m+1
    #
    m = 0 
    while m <= Px-2:   
        S_2[m,0] = Sx[Nuv[m]-1,P1x-1]
        CS_2[m,0] = CSx[Nuv[m]-1,P1x-1]
        m = m+1
    #
    m = 0 
    while m <= Py-2:   
        S_2[Px-1+m,1] = Sy[Nuv[Px-1+m]-1,P2y-1]
        CS_2[Px-1+m,1] = CSy[Nuv[Px-1+m]-1,P2y-1]
        m = m+1 
    #    
    m = 0
    while m <= Px-2:
        S_3[0,m] = Sx[P1x-1,Nuv[m]-1]
        m = m+1
    #
    m = 0
    while m <= Py-2:
        S_3[1,Px-1+m] = Sy[P2y-1,Nuv[Px-1+m]-1]
        m = m+1
    #
    if Py == 1:
        S_422 = Sy[0,0]
        CS_422 = CSy[0,0]
    else:
        S_422 = Sy[P2y-1,P2y-1]
        CS_422 = CSy[P2y-1,P2y-1]    
    #
    S_4 = np.array([[Sx[P1x-1,P1x-1], -1.],[-1., S_422]])
    CS_4 = np.array([[CSx[P1x-1,P1x-1], 0.],[0., CS_422]])       
    #
    # inverse S_4 least square solution
    S_4inv,res,ran,s = np.linalg.lstsq(S_4,np.eye(2), rcond = -1)
    # 
    M_r = np.dot(-S_2,S_4inv)
    M_rtc = np.matrix.getH(M_r)
    #
    CSM_r = np.dot(CS_2,M_rtc)
    CSM_rtc = np.matrix.getH(CSM_r)
    # result S parameters
    rs = np.dot(S_4inv,S_3)
    rs = np.dot(S_2,rs)
    rs = np.add(S_1,-rs)
    # result CS parameters
    rcs = np.dot(CS_4,M_rtc)
    rcs = np.dot(M_r,rcs)
    rcs = np.add(CSM_rtc,rcs)
    rcs = np.add(CSM_r,rcs)
    rcs = np.add(CS_1,rcs)
    #
    return [rs, rcs]
#
##########################################################################
#
def Nnpc_x(Sm,CSm,Px,Py,dsp):
    #
    # port connection at m-port, signal and noise
    #
    #  Sm -> S-matrix m-port
    # CSm -> CS-matrix m-port
    #  Px -> port number from x-port for connection with Py
    #  Py -> port number from x-port for connection with Px
    #        hint: Px > Py !!!
    # dsp -> display paramter to show the port numbers of connection
    #        before and after
    #
    Pm = len(Sm)
    #
    S_1 = np.zeros((Pm-2,Pm-2), dtype = complex)
    CS_1 = np.zeros((Pm-2,Pm-2), dtype = complex)
    #
    S_2 = np.zeros((Pm-2,2), dtype = complex);
    CS_2 = np.zeros((Pm-2,2), dtype = complex);
    #
    S_3 = np.zeros((2,Pm-2), dtype = complex);
    #
    Nx = np.array(list(range(1, Pm+1)))
    Nx0 = Nx
    #
    Nx = np.delete(Nx, Px-1)
    #
    Nx = np.delete(Nx, Py-2)
    #
    Pc = [Px, Py]
    #
    Nuv = np.append(Nx, Pc)
    #
    if dsp == 1:
       N_Tor = len(Nuv)-2
       N_Tor_n = np.array(list(range(1,N_Tor+1)))
       print(' ')
       print('======= connection x port x1 < x2 =============')
       print('== port number x1: ', Px,' ; port number x2: ', Py,' ==')
       print('port numbers start: ', Nx0)
       print('port numbers  open: ', Nx)
       print('port numbers   new: ', N_Tor_n)
       print('===============================================')
    #
    m = 0
    while m <= Pm-3:
        n = 0
        while n <= Pm-3:
            S_1[m,n] = Sm[Nuv[m]-1,Nuv[n]-1]
            CS_1[m,n] = CSm[Nuv[m]-1,Nuv[n]-1]
            n = n+1
        m = m+1
    #
    m = 0 
    while m <= Pm-3:
        S_2[m,0] =  Sm[Nuv[m]-1,Px-1];
        CS_2[m,0] = CSm[Nuv[m]-1,Px-1];
        S_2[m,1] =  Sm[Nuv[m]-1,Py-1];
        CS_2[m,1] = CSm[Nuv[m]-1,Py-1]; 
        m = m+1
    #
    m = 0 
    while m <= Pm-3:
        S_3[0,m] = Sm[Px-1,Nuv[m]-1];
        S_3[1,m] = Sm[Py-1,Nuv[m]-1];
        m = m+1
    #
    S_4 = [[Sm[Px-1,Px-1], Sm[Px-1,Py-1]-1.],
           [Sm[Py-1,Px-1]-1., Sm[Py-1,Py-1]]];
    #
    CS_4 = [[CSm[Px-1,Px-1], CSm[Px-1,Py-1]],
            [CSm[Py-1,Px-1], CSm[Py-1,Py-1]]];
    #
    # inverse S_4 least square solution
    S_4inv,res,ran,s = np.linalg.lstsq(S_4,np.eye(2), rcond = -1)
    # 
    M_r = np.dot(-S_2,S_4inv)
    M_rtc = np.matrix.getH(M_r)
    #
    CSM_r = np.dot(CS_2,M_rtc)
    CSM_rtc = np.matrix.getH(CSM_r)
    # Ergebnis S parameters
    rs = np.dot(S_4inv,S_3)
    rs = np.dot(S_2,rs)
    rs = np.add(S_1,-rs)
    # ergebnis CS parameters
    rcs = np.dot(CS_4,M_rtc)
    rcs = np.dot(M_r,rcs)
    rcs = np.add(CSM_rtc,rcs)
    rcs = np.add(CSM_r,rcs)
    rcs = np.add(CS_1,rcs)
    #
    return [rs, rcs]