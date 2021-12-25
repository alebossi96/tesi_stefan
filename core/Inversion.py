# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 20:12:45 2020

@author: ŠUŠNJAR
"""
import numpy as np
import core.SignalsAndSystems as sis
import subprocess

N_UNKNOWNS = 6

class Inversion(object):
    def __init__(self, W1, W2, signal):
        self.W1 = W1
        self.W2 = W2
        self.signal = signal
    
    def reconstruct_spectra_lsm(self): #When sensitivity matrix is W(x)
        Wt = np.matrix([self.W1, self.W2])
        W = Wt.transpose()
        S1 = sis.Spectrum("First layer spectrum")
        S2 = sis.Spectrum("Second layer spectrum")
        S1.n_points = self.signal.n_points
        S2.n_points = self.signal.n_points
        S1.l = np.zeros(S1.n_points,dtype=np.float)
        S2.l = np.zeros(S2.n_points,dtype=np.float)
        S1.s = np.zeros(S1.n_points,dtype=np.float)
        S2.s = np.zeros(S2.n_points,dtype=np.float)
        for l in range(0,self.signal.n_points):
            S1.l[l] = self.signal.l[l]
            S2.l[l] = self.signal.l[l]
            A = Wt*W
            b = Wt*np.matrix(self.signal.s[l]).transpose()
            x = np.linalg.solve(A, b)
            S1.s[l] = x[0]
            S2.s[l] = x[1]
        return S1, S2
    
    def reconstruct_spectra_lsm_wavelength_dependent(self): #When sensitivity matrix is W(lambda,x)
        S1 = sis.Spectrum("First layer spectrum")
        S2 = sis.Spectrum("Second layer spectrum")
        S1.n_points = self.signal.n_points
        S2.n_points = self.signal.n_points
        S1.l = np.zeros(S1.n_points,dtype=np.float)
        S2.l = np.zeros(S2.n_points,dtype=np.float)
        S1.s = np.zeros(S1.n_points,dtype=np.float)
        S2.s = np.zeros(S2.n_points,dtype=np.float)
        for l in range(0,self.signal.n_points):
            Wt = np.matrix([self.W1[l], self.W2[l]])
            W = Wt.transpose()
            S1.l[l] = self.signal.l[l]
            S2.l[l] = self.signal.l[l]
            A = Wt*W
            b = Wt*np.matrix(self.signal.s[l]).transpose()
            x = np.linalg.solve(A, b)
            S1.s[l] = x[0]
            S2.s[l] = x[1]
        return S1, S2

def norm_square(n, v):
    """
    

    Parameters
    ----------
    n : INT
        LENGTH OF VECTOR.
    v : COLUMN VECTOR (MATRIX WITH DIMENSIONS nx1) OF FLOAT
        VECTOR WHICH NORM IS TO BE CALCULATED.

    Returns
    -------
    s : FLOAT
        SQUARE OF THE VECTOR 2-NORM OF VECTOR v.

    """
    s = 0
    for i in range(0,n):
        s+=v[i,0]**2
    return s

def iterative_method(n_gates, tmin, tmax, di, b, x0):
   
   N_ITER_MAX = 100
   it = 0
   #x0 = [0.005, 0.005, 1., 1., 0.001, 0.001]
   xk = np.matrix(x0).transpose() #Initial attempt
   lk = 1e-6 #Initial 'dumping' parameter lambda(k), for k=0
   #nu = 1.5 #Initial coefficient \nu (should be greater than 1)
   arguments_list_base = "TwoLayerJacobian.exe "+str(di.z0)+" "+str(di.z1)+" "+str(di.rec[0])+" "+str(di.R)+" "+str(di.n0)+" "+str(di.n1)+" "+str(di.ne)+" "+str(di.c)+" "+str(di.dt)+" "+str(di.tmin)+" "+str(di.dr)+" "+str(di.rmin)+" "+str(di.n_tpsf)+" "+str(n_gates)+" "+str(tmin)+" "+str(tmax)+" "
   arguments_list = arguments_list_base+str(x0[0])+" "+str(x0[1])+" "+str(x0[2])+" "+str(x0[3])+" "+str(x0[4])+" "+str(x0[5])
   subprocess.call(arguments_list,shell=True)
   
   fr = open("Jacobian.txt","r")
   J = np.matrix(np.zeros((n_gates,N_UNKNOWNS)),np.float)
   f = np.matrix(np.zeros((n_gates,1)),np.float)
   for r in range(0,n_gates):
       line = fr.readline()
       for c in range(0,N_UNKNOWNS):
           ind = line.find(' ')
           J[r,c] = float(line[:ind])
           line = line[ind+1:]
       f[r,0] = float(line[:-1])
   fr.close()
   dxk = np.linalg.solve(J.transpose()*J+lk*np.matrix(np.identity(N_UNKNOWNS,np.float)),J.transpose()*(b-f))
   xk = xk+dxk
   it+=1
   # print(xk)
   # print("\n")
   # print(dxk)
   # print("\n")
   fik = norm_square(N_UNKNOWNS, b-f)
   rel_err = norm_square(N_UNKNOWNS, b-f)/norm_square(N_UNKNOWNS,b)

   while(it<N_ITER_MAX):
       print("Iteration step:%d %e %e\n" % (it, lk, fik))
       # print(xk)
       # print("\n")
       # print(dxk)
       # print("\n")
       fr = open("Jacobian.txt","r")
       #J = np.matrix(np.zeros((n_gates,N_UNKNOWNS)),np.float)
       #f = np.matrix(np.zeros((n_gates,1)),np.float)
       for r in range(0,n_gates):
           line = fr.readline()
           for c in range(0,N_UNKNOWNS):
               ind = line.find(' ')
               J[r,c] = float(line[:ind])
               line = line[ind+1:]
           f[r,0] = float(line[:-1])
       fr.close()
       dxk = np.linalg.solve(J.transpose()*J+lk*np.matrix(np.identity(N_UNKNOWNS,np.float)),J.transpose()*(b-f))
       fik_prev = fik
       fik = norm_square(n_gates, b-f)
       
       if (fik<=fik_prev):
           lk = 0.8*lk
           xk = xk+dxk
           it+=1
           arguments_list = arguments_list_base+str(xk[0,0])+" "+str(xk[1,0])+" "+str(xk[2,0])+" "+str(xk[3,0])+" "+str(xk[4,0])+" "+str(xk[5,0])
           subprocess.call(arguments_list,shell=True)
           rel_err = fik/norm_square(n_gates,b)
       else:
           lk = 2*lk
       if (rel_err<1e-6):
           break
   print("\nIteration process finished.\n")
   return xk
           
def iterative_method_Broyden(n_gates, tmin, tmax, di, b, x0):
    N_ITER_MAX = 100
    it = 0
    xk = np.matrix(x0).transpose() #Initial attempt
    arguments_list_base = "TwoLayerJacobian.exe "+str(di.z0)+" "+str(di.z1)+" "+str(di.rec[0])+" "+str(di.R)+" "+str(di.n0)+" "+str(di.n1)+" "+str(di.ne)+" "+str(di.c)+" "+str(di.dt)+" "+str(di.tmin)+" "+str(di.dr)+" "+str(di.rmin)+" "+str(di.n_tpsf)+" "+str(n_gates)+" "+str(tmin)+" "+str(tmax)+" "
    arguments_list = arguments_list_base+str(x0[0])+" "+str(x0[1])+" "+str(x0[2])+" "+str(x0[3])+" "+str(x0[4])+" "+str(x0[5])
    subprocess.call(arguments_list,shell=True)
    
    fr = open("Jacobian.txt","r")
    J = np.matrix(np.zeros((n_gates,N_UNKNOWNS)),np.float)
    f = np.matrix(np.zeros((n_gates,1)),np.float)
    for r in range(0,n_gates):
        line = fr.readline()
        for c in range(0,N_UNKNOWNS):
            ind = line.find(' ')
            J[r,c] = float(line[:ind])
            line = line[ind+1:]
        f[r,0] = float(line[:-1])
    fr.close()
    Qk = J #still k=0
    while(it<N_ITER_MAX):
        dxk = np.linalg.lstsq(Qk,b-f,rcond=None)
        dxk = dxk[0]
        xk = xk+dxk
        arguments_list = arguments_list_base+str(xk[0,0])+" "+str(xk[1,0])+" "+str(xk[2,0])+" "+str(xk[3,0])+" "+str(xk[4,0])+" "+str(xk[5,0])
        subprocess.call(arguments_list,shell=True)
        fr = open("Jacobian.txt","r")
        f_prev = f
        for r in range(0,n_gates):
            line = fr.readline()
            for c in range(0,N_UNKNOWNS):
                ind = line.find(' ')
                J[r,c] = float(line[:ind])
                line = line[ind+1:]
            f[r,0] = float(line[:-1])
        fr.close()
        dfk = f - f_prev
        Qk = Qk + ((dfk-Qk*dxk)*(dxk.transpose()))/(dxk.transpose()*dxk)
        it+=1
        rel_err = norm_square(n_gates,b-f)/norm_square(n_gates, b)
        print("Iteration step %d: rel_err=%e\n" % (it,rel_err))
        if (rel_err<1e-4):
            break
    return xk

def guess_nearest_solution(sig):
    N_A_2 = 50 #2*N_A_2+1 is total number of points for ua2 guessing
    N_S_2 = 50 #2*N_S_2+1 is total number of points for ud2 guessing
    n_gates = 10 #All these could be sent as parameters when calling subprocess
    rows = (2*N_A_2+1)*(2*N_S_2+1)
    Wg1 = np.zeros((rows,n_gates),dtype=np.float)
    Wg2 = np.zeros((rows,n_gates),dtype=np.float)
    ua2 = np.zeros(rows,dtype=np.float)
    us2 = np.zeros(rows,dtype=np.float)
    S1 = sis.Spectrum("First layer spectrum")
    S2 = sis.Spectrum("Second layer spectrum")
    S1.n_points = sig.n_points
    S2.n_points = sig.n_points
    S1.l = np.zeros(S1.n_points,dtype=np.float)
    S2.l = np.zeros(S2.n_points,dtype=np.float)
    S1.s = np.zeros(S1.n_points,dtype=np.float)
    S2.s = np.zeros(S2.n_points,dtype=np.float) 
    ua2_spectrum = sis.Spectrum("Second layer absorption") #Corresponding
    us2_spectrum = sis.Spectrum("Second layer scattering") #spectra
    ua2_spectrum.n_points = sig.n_points
    us2_spectrum.n_points = sig.n_points
    ua2_spectrum.l = np.zeros(ua2_spectrum.n_points,dtype=np.float)
    us2_spectrum.l = np.zeros(us2_spectrum.n_points,dtype=np.float)
    ua2_spectrum.s = np.zeros(ua2_spectrum.n_points,dtype=np.float)
    us2_spectrum.s = np.zeros(us2_spectrum.n_points,dtype=np.float)
    fr = open("Sensitivity.txt","r")
    print("Reading from file started.\n")
    for i in range(0,rows):
        line = fr.readline()
        ind = line.find(' ')
        ua2[i] = float(line[:ind])
        line = line[ind+1:]
        ind = line.find(' ')
        us2[i] = float(line[:ind])
        line = line[ind+1:]
        for j in range(0,n_gates-1): #Reading sensitivity matrix from file
            ind = line.find(' ')
            Wg1[i][j] = float(line[:ind])
            line = line[ind+1:]
            ind = line.find(' ')
            Wg2[i][j] = float(line[:ind])
            line = line[ind+1:]
        ind = line.find(' ')
        Wg1[i][j] = float(line[:ind])
        line = line[ind+1:]
        Wg2[i][j] = float(line[:-1])
    fr.close() #End of reading
    print("Reading from file finished.\n")
    #Solving at each wavelength
    for l in range(0, sig.n_points):
        S1.l[l] = sig.l[l]
        S2.l[l] = sig.l[l]
        ua2_spectrum.l[l] = sig.l[l]
        us2_spectrum.l[l] = sig.l[l]
        norm_min = 999999 #Minimum distance from real signal
        print("Wavelength: %d (%lf)   " % (l+1,S1.l[l]))
        for i in range(0,rows):
            Wt = np.matrix([Wg1[i], Wg2[i]])
            W = Wt.transpose()
            A = Wt*W
            b = Wt*np.matrix(sig.s[l]).transpose()
            x = np.linalg.solve(A, b)
            norm_current = 0
            for j in range(0, n_gates):
                norm_current += (sig.s[l][j]-Wg1[i][j]*x[0]-Wg2[i][j]*x[1])**2
            #print("Err=%e\n" % norm_current)
            if (norm_current<norm_min):
                ua2_spectrum.s[l] = ua2[i]
                us2_spectrum.s[l] = us2[i]
                S1.s[l] = x[0]
                S2.s[l] = x[1]
                norm_min = norm_current
    #End of for loop for solving at each wavelength
    return ua2_spectrum, us2_spectrum, S1, S2           
            
            
        
    

   
    
   
                
        
