# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:39:32 2020

@author: ŠUŠNJAR
"""
from math import*
import numpy as np

def read_data_line_from_file(file, key):
    """

    Parameters
    ----------
    file : FILE
        FILE POINTER, ASSOCIATED WITH OPENED INPUT FILE IN THE READ MODE ('R').
    key : STRING
        STRING THAT DESCRIBES THE VARIABLE BEING READ FROM THE INPUT FILE.

    Returns
    -------
    (0, string) IF READING WAS SUCCESFUL OR (-1, string) IF THE FORMAT WAS NOT EXPECTED.

    """														  
    string = file.readline()
    ind = string.find('=')
    if (ind > 0 and key == string[:ind]):
        string = string[ind+1:-1]
        return 0,string
    else:
        return -1,string


def dAA_cyl (an12):
    A=0.0
    if (an12 > 1.0):
        A=504.332889-2641.00214*an12+5923.699064*an12*an12+-7376.355814*an12*an12*an12+5507.5304*pow(an12,4)-2463.357945*pow(an12,5.)+610.956547*pow(an12,6.)-64.8047*pow(an12,7.)
    elif (an12 == 1.):
        A=1.0
    else: 
        A=3.084635-6.531194*an12+8.357854*an12*an12-5.082751*an12*an12*an12+1.171382*pow(an12,4.)
    return A

def ftd_cyl_2_r_r(square, kz0, f, df, di):
    # Funzione che calcola la funzione e la derivata dell'equazione trascendente
    # Caso considerato = Reale - Reale
    # Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
    # Version modified for Python
          
    kz1=sqrt((pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.)))
    dkz1=kz0/kz1*(di.ud1/di.ud0)*(di.n1/di.n0)
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)
    ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*d1*kz1)
    ga0=(2*di.A0e*d0*kz0)
    dga0=2*di.A0e*d0
    dga1=-dkz1*(di.z0+di.z1)-2*di.A1e*d1*dkz1
    f=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-(3*di.ud1/kz1)*tan(kz1*di.z0+ga1)*pow((di.n0/di.n1),2.)
    df=-(3*di.ud0/pow(kz0,2.))*tan(kz0*di.z0+ga0)+(3*di.ud1/pow(kz1,2.))*tan(kz1*di.z0+ga1)*dkz1*pow((di.n0/di.n1),2.)+(3*di.ud0/kz0)/pow((cos(kz0*di.z0+ga0)),2.)*(di.z0+dga0)-(3*di.ud1/kz1)/pow((cos(kz1*di.z0+ga1)),2.)*(dkz1*di.z0+dga1)*pow((di.n0/di.n1),2.)
    return f, df

def ftd_cyl_2_r_i(square, kz0, f, df, di):
    # Funzione che calcola la funzione e la derivata dell'equazione trascendente
    # Caso considerato = Reale - Immaginario
    # Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
    # Version modified for Python
    #     kz1 is immaginary
    #     f and df are the values of the function and of his derivative with respect kz0
    
    kz1_2=(pow(kz0,2)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.))
    mod_kz1=sqrt(abs(-kz1_2))
    dmod_kz1=-kz0/mod_kz1*(di.ud1/di.ud0)*(di.n1/di.n0)
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)      
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1)
    ga0=(2*di.A0e*d0*kz0)
    dga0=2*di.A0e*d0
    dga1=-dmod_kz1*(di.z0+di.z1)-2*di.A1e*d1*dmod_kz1
    arg=mod_kz1*di.z0+ga1
    darg=dmod_kz1*di.z0+dga1	   
    f=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-(3*di.ud1/mod_kz1)*tanh(arg)*pow((di.n0/di.n1),2.)
    df=-(3*di.ud0/pow(kz0,2.))*tan(kz0*di.z0+ga0)+(3*di.ud1/pow(mod_kz1,2.))*dmod_kz1*tanh(arg)*pow((di.n0/di.n1),2.)+(3*di.ud0/kz0)/pow((cos(kz0*di.z0+ga0)),2.)*(di.z0+dga0)-(3*di.ud1/mod_kz1)/pow(cosh(arg),2.)*darg*pow((di.n0/di.n1),2.)
    return f, df

def ftd_cyl_2_i_r(square, kz0, f, df, di):
    # Funzione che calcola la funzione e la derivata dell'equazione trascendente
    # Caso considerato = Reale - Immaginario
    # Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
    # Version modified for Python
    #     kz0 is immaginary
    #     La variabile kz0 rappresenta in questo caso il modulo |kn0| della variabile kn0
    #     f and df are the values of the function and of his derivative with respect |kn0| 
    kz1_2=(-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.))
    mod_kz1=sqrt(abs(kz1_2))
    dmod_kz1=-kz0/mod_kz1*(di.ud1/di.ud0)*(di.n1/di.n0)
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)      
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1)
    ga0=(2*di.A0e*d0*kz0)
    dga0=2*di.A0e*d0
    dga1=-dmod_kz1*(di.z0+di.z1)-2*di.A1e*d1*dmod_kz1
    arg=mod_kz1*di.z0+ga1
    darg=dmod_kz1*di.z0+dga1	   
    f=(3*di.ud0/kz0)*tanh(kz0*di.z0+ga0)-(3*di.ud1/mod_kz1)*tan(arg)*pow((di.n0/di.n1),2.)
    df=-(3*di.ud0/pow(kz0,2.))*tanh(kz0*di.z0+ga0)+(3*di.ud1/pow(mod_kz1,2.))*dmod_kz1*tan(arg)*pow((di.n0/di.n1),2.)+(3*di.ud0/kz0)/pow((cosh(kz0*di.z0+ga0)),2.)*(di.z0+dga0)-(3*di.ud1/mod_kz1)/pow(cos(arg),2.)*darg*pow((di.n0/di.n1),2.)  
    return f, df

def ft_cyl_2_r_i(square, kz0, di):
# Funzione che calcola la funzione e la derivata dell'equazione trascendente
# Caso considerato = Reale - Immaginario
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
# Version modified for Python
#     kz1 is immaginary
#     f and df are the value of the function with respect kz0
    kz1_2=(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.))
    mod_kz1=sqrt(abs(-kz1_2))   #   valore assoluto inserito per precauzione
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)       
    ga0=(2*di.A0e*d0*kz0)
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1)
    arg=mod_kz1*di.z0+ga1   
    f_out=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-(3*di.ud1/mod_kz1)*tanh(arg)*pow((di.n0/di.n1),2.)
    return f_out

def ft_cyl_2_r_r(square, kz0, di):
# Funzione che calcola la funzione dell'equazione trascendente
# Caso considerato = Reale - Reale
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab      
# Version modified for Python
    kz1=sqrt(abs(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.)))
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)      
    ga0=(2*di.A0e*d0*kz0)
    ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*d1*kz1)
    f_out=(3*di.ud0/kz0)*tan(kz0*di.z0+ga0)-(3*di.ud1/kz1)*tan(kz1*di.z0+ga1)*pow((di.n0/di.n1),2.)
    return f_out

def ft_cyl_2_i_r(square, kz0, di):
# Funzione che calcola la funzione e la derivata dell'equazione trascendente
# Caso considerato = Reale - Immaginario
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
# Version modified for Python
#     kz0 is immaginary
#     La variabile kz0 rappresenta in questo caso il modulo |kn0| della variabile kn0
#     f and df are the values of the function and of his derivative with respect |kn0|
    kz1_2=(-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(square)*(di.ud1/di.ud0*di.n1/di.n0-1.))
    mod_kz1=sqrt(abs(kz1_2))   #   valore assoluto dentro la radice inserito per precauzione
    d0=1./(3*di.ud0)            # per problemi di precisione si possono avere valori piccolissi e negativi di kz1_2
    d1=1./(3*di.ud1)      
    ga0=(2*di.A0e*d0*kz0)
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1)
    arg=mod_kz1*di.z0+ga1   
    f_out=(3*di.ud0/kz0)*tanh(kz0*di.z0+ga0)-(3*di.ud1/mod_kz1)*tan(arg)*pow((di.n0/di.n1),2.)
    return f_out

def r_cyl_2_r_r(t, kz0, kj, di):      
# Funzione che calcola la riflettanza da un mezzo a due strati
# Caso considerato = Reale - Reale
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
# Version modified for Python
    j1=bessj1(di.R*kj)
    j0=bessj0(di.ro*kj)
    # pi=2.d0*dasin(1.d0)
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)      
    kz1=sqrt(abs(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(pow(kj,2.))*(di.ud1/di.ud0*di.n1/di.n0-1.)))
    ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*d1*kz1)
    ga0=(2*di.A0e*d0*kz0)       
    arg1=sin(kz1*di.z0+ga1)   
    b1=sin(kz0*di.z0+ga0)/arg1*pow((di.n1/di.n0),2.)
    Nl=di.z0/2.+ga0/(2.*kz0)-pow(b1,2.)/(2*kz1)*(kz1*di.z0+ga1)*(di.n0/di.n1)-(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))+(pow(b1,2.)/(4*kz1))*(sin(2*(kz1*di.z0+ga1)))*(di.n0/di.n1)
    # Il contributo del termine radiale // e' inserito nella fil2_cyl
    f_out=d0*di.v0*kz0*cos(ga0)*sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua0*di.v0*t)*exp(-(pow(kj,2.)+pow(kz0,2.))*d0*di.v0*t)
    return f_out

def r_cyl_2_r_i(t, kz0, kj, di):
# Funzione che calcola la riflettanza da un mezzo a due strati
# Caso considerato = Reale - Immaginario
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab      
# Version modified for Python
    j1=bessj1(di.R*kj);
    j0=bessj0(di.ro*kj);
    # pi=2.d0*dasin(1.d0)
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)
    kz1_2=((pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(pow(kj,2.))*(di.ud1/di.ud0*di.n1/di.n0-1.)))
    mod_kz1=sqrt(-kz1_2)
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1)
    ga0=(2*di.A0e*d0*kz0)       
    arg1=(mod_kz1*di.z0+ga1);
    #	  if (abs(2*arg1) > 300.):
    #       b1=d0
      
    b1=sin(kz0*di.z0+ga0)/sinh(arg1)*pow((di.n1/di.n0),2.)
#      b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
    Nl=di.z0/2.+ga0/(2.*kz0)+pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1)-(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))-(pow(b1,2.)/(4*mod_kz1))*(sinh(2*arg1))*(di.n0/di.n1)
    # Il contributo del termine radiale // e' inserito direttamente nella fil2_r_i
#     +-(pow(b1,2.)/(4*mod_kz1))*((exp(-2*arg1)-exp(2*arg1))/2)
    f_out=d0*di.v0*kz0*cos(ga0)*sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua0*di.v0*t)*exp(-(pow(kj,2.)+pow(kz0,2.))*d0*di.v0*t)
    return f_out

def r_cyl_2_i_r(t, kz0, kj, di):      
# Funzione che calcola la riflettanza da un mezzo a due strati
# Caso considerato = Reale - Immaginario
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
# Version modified for Python
    j1=bessj1(di.R*kj)
    j0=bessj0(di.ro*kj)
    # pi=2.d0*dasin(1.d0)
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)      
    kz1_2=((-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(pow(kj,2.))*(di.ud1/di.ud0*di.n1/di.n0-1.)))
    mod_kz1=sqrt(abs(kz1_2))
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*d1*mod_kz1);
    ga0=(2*di.A0e*d0*kz0);       
    arg1=(mod_kz1*di.z0+ga1);   
    b1=sinh(kz0*di.z0+ga0)/sin(arg1)*pow((di.n1/di.n0),2.);
#  b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
    Nl=-di.z0/2.-ga0/(2.*kz0)-pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1)+(1./(4.*kz0))*(sinh(2*(kz0*di.z0+ga0)))+(pow(b1,2.)/(4*mod_kz1))*(sin(2*arg1))*(di.n0/di.n1) 
   # IL contributo del termine radiale e' inserito direttamente nella fil2_i_r
#     +-(b1**2/(4*mod_kz1))*((dexp(-2*arg1)-dexp(2*arg1))/2)
    f_out=d0*di.v0*kz0*cosh(ga0)*sinh(kz0*3*d0+ga0)/(di.pi*pow(di.R,2)*pow(j1,2.))*j0/Nl*exp(-di.ua0*di.v0*t)*exp((-pow(kj,2.)+pow(kz0,2.))*d0*di.v0*t);   
    return f_out

def t_cyl_2_r_r(t, kz0, kj, di):
# Funzione che calcola la trasmittanza attraverso un mezzo a due strati
# Caso considerato = Reale - Reale
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
# Version modified for Python
    j1=bessj1(di.R*kj)
    j0=bessj0(di.ro*kj)
    # pi=2.d0*dasin(1.d0)
    d0=1./(3*di.ud0)
    d1=1./(3*di.ud1)
    kz1=sqrt(abs(pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(pow(kj,2.))*(di.ud1/di.ud0*(di.n1/di.n0)-1.)))
    ga1=-kz1*(di.z0+di.z1)-(2*di.A1e*(1./(3*di.ud1))*kz1)
    ga0=(2*di.A0e*d0*kz0)
    arg1=sin(kz1*di.z0+ga1)
    b1=sin(kz0*di.z0+ga0)/arg1*pow((di.n1/di.n0),2.)
    Nl=di.z0/2.+ga0/(2.*kz0)-pow(b1,2.)/(2*kz1)*(kz1*di.z0+ga1)*(di.n0/di.n1)-(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))+(pow(b1,2.)/(4*kz1))*(sin(2*(kz1*di.z0+ga1)))*(di.n0/di.n1)
    # Il contributo del termine radiale e' inserito nella fil2_cyl
    f_out=-d1*di.v1*di.n1/di.n0*b1*kz1*cos(kz1*(di.z0+di.z1)+ga1)*sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua1*di.v1*t)*exp(-(pow(kj,2.)+pow(kz1,2.))*d1*di.v1*t)
    return f_out

def t_cyl_2_r_i(t, kz0, kj, di):      
# Funzione che calcola la trasmittanza attraverso un mezzo a due strati
# Caso considerato = Reale - Immaginario
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
# Version modified for Python
    j1=bessj1(di.R*kj);
    j0=bessj0(di.ro*kj);
    # pi=2.d0*dasin(1.d0)
    d0=1./(3*di.ud0);
    d1=1./(3*di.ud1);      
    kz1_2=((pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(pow(kj,2.))*(di.ud1/di.ud0*(di.n1/di.n0)-1.)))
    mod_kz1=sqrt(-kz1_2)
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*(1./(3*di.ud1))*mod_kz1)
    ga0=(2*di.A0e*d0*kz0)       
    arg1=(mod_kz1*di.z0+ga1)   
    b1=sin(kz0*di.z0+ga0)/sinh(arg1)*pow((di.n1/di.n0),2.)
    # b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
    Nl=di.z0/2.+ga0/(2.*kz0)+pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1)-(1./(4.*kz0))*(sin(2*(kz0*di.z0+ga0)))-(pow(b1,2.)/(4*mod_kz1))*(sinh(2*arg1))*(di.n0/di.n1)
    # Il contributo del termine radiale // e' inserito direttamente nella fil2_r_i
    #     +-(b1**2/(4*mod_kz1))*((dexp(-2*arg1)-dexp(2*arg1))/2)
    f_out=-d1*di.v1*di.n1/di.n0*mod_kz1*b1*cosh(mod_kz1*(di.z0+di.z1)+ga1)*sin(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua1*di.v1*t)*exp(-(pow(kj,2.)-pow(mod_kz1,2.))*d1*di.v1*t)
    return f_out

def t_cyl_2_i_r (t, kz0, kj, di):      
# Funzione che calcola la trasmittanza attraverso un mezzo a due strati
# Caso considerato = Immaginario - Reale
# Versione in C scritta a partire dal programma Fortran Cyl_lay2_turbo_fab
# Version modified for Python
    j1=bessj1(di.R*kj);
    j0=bessj0(di.ro*kj);
    # pi=2.d0*dasin(1.d0)
    d0=1./(3*di.ud0);
    d1=1./(3*di.ud1);    
    kz1_2=((-pow(kz0,2.)*(di.ud1/di.ud0)*(di.n1/di.n0)+(3*di.ud1)*(di.ua0*(di.n1/di.n0)-di.ua1)+(pow(kj,2.))*(di.ud1/di.ud0*(di.n1/di.n0)-1.)))
    mod_kz1=sqrt(abs(kz1_2))
    ga1=-mod_kz1*(di.z0+di.z1)-(2*di.A1e*(1./(3*di.ud1))*mod_kz1)
    ga0=(2*di.A0e*d0*kz0)
    arg1=(mod_kz1*di.z0+ga1)
    b1=sinh(kz0*di.z0+ga0)/sin(arg1)*pow((di.n1/di.n0),2.)
    # b1=2*dsin(kz0*z0+ga0)/(dexp(-arg1)-dexp(arg1))
    Nl=-di.z0/2.-ga0/(2.*kz0)-pow(b1,2.)/(2*mod_kz1)*(arg1)*(di.n0/di.n1)+(1./(4.*kz0))*(sinh(2*(kz0*di.z0+ga0)))+(pow(b1,2.)/(4*mod_kz1))*(sin(2*arg1))*(di.n0/di.n1)
    # IL contributo del termine radiale e' inserito direttamente nella fil2_i_r
    #    +-(b1**2/(4*mod_kz1))*((dexp(-2*arg1)-dexp(2*arg1))/2)
    f_out=-d1*di.v1*di.n1/di.n0*mod_kz1*b1*cos(mod_kz1*(di.z0+di.z1)+ga1)*sinh(kz0*3*d0+ga0)/(di.pi*pow(di.R,2.)*pow(j1,2.))*j0/Nl*exp(-di.ua1*di.v1*t)*exp((-pow(kj,2.)-pow(mod_kz1,2.))*d1*di.v1*t)
    return f_out

def bessj0(x):
    ax = abs(x)
    if (ax < 8.0):
        y=x*x;
        ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))))
        ans2=57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))))
        ans=ans1/ans2
    else:
        z=8.0/ax
        y=z*z
        xx=ax-0.785398164
        ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
        ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934935152e-7)))
        ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    return ans


def bessj1(x):
    ax = abs(x)
    if (ax < 8.0):
        y=x*x
        ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))))
        ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))))
        ans=ans1/ans2
    else:
        z=8.0/ax
        y=z*z
        xx=ax-2.356194491
        ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))))
        ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)))
        ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2)
        if (x < 0.0): 
            ans = -ans
    return ans


MAXIT = 3000
def rtsafe_cyl_2(funcd, x1, x2, xacc, square, di):
    fl = 0
    fh = 0
    df = 0
    f = 0
    (fl,df)=funcd(square,x1,fl,df,di)
    (fh,df)=funcd(square,x2,fh,df,di)
    if ((fl > 0.0 and fh > 0.0) or (fl < 0.0 and fh < 0.0)):
#		nrerror("Root must be bracketed in rtsafe");
        pass
    if (fl == 0.0):
        return x1
    if (fh == 0.0):
        return x2
    if (fl < 0.0):
        xl=x1
        xh=x2
    else:
        xh=x1
        xl=x2
        
    rts=0.5*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    (f,df)=funcd(square,rts,f,df,di)
    for j in range(1,MAXIT+1):
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) or (abs(2.0*f) > abs(dxold*df))):
            dxold=dx
            dx=0.5*(xh-xl)
            rts=xl+dx
            if (xl == rts):
                return rts
        else:
            dxold=dx
            dx=f/df
            temp=rts
            rts -= dx
            if (temp == rts):
                return rts;
            if (abs(dx) < xacc):
                return rts
        (f,df)=funcd(square,rts,f,df,di)
        if (f < 0.0):
            xl=rts;
        else:
            xh=rts
#	nrerror("Maximum number of iterations exceeded in rtsafe");
    return 0.0


def zbrak_cyl_2(fx, x1, x2, n, xb1, xb2, nb, di, square):
    nbb=0;
    dx=(x2-x1)/n;
    x=x1
    fp=fx(square,x,di)
    for i in range(1,n+1):
        x = x+dx
        fc=fx(square,x,di)
#	fc=(*fx)(square,x =x1+dx*i,di); 
        if (fc*fp <= 0.0):
            nbb=nbb+1
            xb1[nbb-1]=x-dx
            xb2[nbb-1]=x
            if(nb == nbb): 
                return nb
        fp=fc;
    nb = nbb;
    return nb
        
def f_roots(di, d_o):
    #  Funzione che calcola le radici dell'equazione trascendente
    # di - Data input, d_o - Raw data output   
    NBB = 6000
    N_D = 2000
    fx1 = 0.
    dfx1 = 0.
    fx2 = 0.
    dfx2 = 0.
    zb1_r_i = np.zeros(NBB,dtype=np.float)
    zb1_i_r = np.zeros(NBB,dtype=np.float)
    zb2_r_i = np.zeros(NBB,dtype=np.float)
    zb2_i_r = np.zeros(NBB,dtype=np.float)
    xb1 = np.zeros(N_D,dtype=np.float)
    xb2 = np.zeros(N_D,dtype=np.float)
    discontinuita_1 = np.zeros(N_D,dtype=np.float)
    discontinuita_2 = np.zeros(N_D,dtype=np.float)
    idisc_t=0
    pi_2=di.pi/2.
    d0_0=1./(3*di.ud0)
    d1_1=1./(3*di.ud1)   
    

    for ii in range(0,di.n_kj):
        square=pow(d_o.kappa_j[ii],2.);  # kappa_j(ii)**2
        ref=square*(d1_1/di.n1-d0_0/di.n0)+(di.ua1/di.n1-di.ua0/di.n0);  # (kappa_j(ii)**2)*(D1/n1-D0/n0)+ua1/n1-ua0/n0
        #     Rispetto alla teoria del manoscritto il termine ref e' uguale a -C*D1/n1
        #     in the following if control we are ckecking which are the smallest
        #     values of the k coefficients 
        if (ref < 0.):
            kz1min=sqrt(-ref/d1_1*di.n1)
            kz0min=di.precisione
        elif (ref > 0.):
            kz0min=sqrt(ref/d0_0*di.n0)
            kz1min=di.precisione
        else:
            kz0min=di.precisione
            kz1min=di.precisione
        nzz_r_i=0
        nzz_i_r=0 	  
#-----------------------------------------------------------*
#     Bracketing of immaginary roots	
        if (ref >= 0.):
            # kz0 real, kz1 immaginary
            M0=(int)((sqrt(ref/d0_0*di.n0))*(di.z0+2*di.A0e*d0_0)/di.pi)
            if (M0 > 0):
                x1=M0*di.pi/((di.z0+2*di.A0e*d0_0))+di.precisione
                x2=sqrt(ref/d0_0*di.n0)-di.precisione
                nb=6000    
                nb = zbrak_cyl_2(ft_cyl_2_r_i,x1,x2,di.nh,xb1,xb2,nb,di,square)
                if (nb > 0):
                    for j in range(0,nb):
                        nzz_r_i=nzz_r_i+1
                        zb1_r_i[nzz_r_i-1]=xb1[j]
                        zb2_r_i[nzz_r_i-1]=xb2[j]
             
                for i in range(1,M0+1):
                    x1=(2*i-1)*di.pi/(2*(di.z0+2*di.A0e*d0_0))
                    x2=(2*i)*di.pi/(2*(di.z0+2*di.A0e*d0_0))
                    nzz_r_i=nzz_r_i+1
                    zb1_r_i[nzz_r_i-1]=x1+di.precisione
                    zb2_r_i[nzz_r_i-1]=x2-di.precisione
	       
            elif (M0 == 0 ):
                if (kz0min*(di.z0+2*di.A0e*d0_0) > pi_2):
                    x1=(pi_2)/(di.z0+2*di.A0e*d0_0)+di.precisione
                    x2=sqrt(ref/d0_0*di.n0)-di.precisione
                    nb=6000
                    nb = zbrak_cyl_2(ft_cyl_2_r_i,x1,x2,di.nh,xb1,xb2,nb,di,square)
                    if (nb > 0):
                        for j in range(0,nb):
                            nzz_r_i=nzz_r_i+1
                            zb1_r_i[nzz_r_i-1]=xb1[j]
                            zb2_r_i[nzz_r_i-1]=xb2[j]
           
        elif (ref < 0.):    # kz0 immaginary, kz1 real
            M0=(int)((sqrt(-ref/d1_1*di.n1))*(di.z1+2*di.A1e*d1_1)/di.pi)
            if (M0 > 0):
                x1=di.precisione
                x2=sqrt(-pow((M0*di.pi/((di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*(di.n0/di.n1)-ref/d0_0*di.n0)-di.precisione
                nb=6000    
                nb = zbrak_cyl_2(ft_cyl_2_i_r,x1,x2,di.nh,xb1,xb2,nb,di,square)   
                if (nb > 0):
                    for j in range(0,nb):
                        nzz_i_r=nzz_i_r+1 		    
                        zb1_i_r[nzz_i_r-1]=xb1[j]
                        zb2_i_r[nzz_i_r-1]=xb2[j]
          
                for i in range(1,M0+1):
                    x1=sqrt(-pow(((2*i)*di.pi/(2*(di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*di.n0/di.n1-ref/d0_0*di.n0)+di.precisione
                    x2=sqrt(-pow(((2*i-1)*di.pi/(2*(di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*di.n0/di.n1-ref/d0_0*di.n0)-di.precisione
                    nzz_i_r=nzz_i_r+1
                    zb1_i_r[nzz_i_r-1]=x1+di.precisione
                    zb2_i_r[nzz_i_r-1]=x2-di.precisione
            elif(M0 == 0):
                if (kz1min*(di.z1+2*di.A1e*d1_1) > pi_2):
                    x1=di.precisione
                    x2=sqrt(-pow((di.pi/(2*(di.z1+2*di.A1e*d1_1))),2.)*d1_1/d0_0*di.n0/di.n1-ref/d0_0*di.n0)-di.precisione
                    nb=6000    
                    nb = zbrak_cyl_2(ft_cyl_2_i_r,x1,x2,di.nh,xb1,xb2,nb,di,square)
                    if (nb > 0):
                        for j in range(0,nb):
                            nzz_i_r=nzz_i_r+1        
                            zb1_i_r[nzz_i_r-1]=xb1[j]
                            zb2_i_r[nzz_i_r-1]=xb2[j]
#-----------------------------------------------------------*   
#-----------------------------------------------------------
#     Calcolo delle radici immaginarie Caso ref>0 ! kz0 real, kz1 immaginary
        if (ref > 0.):
            i_r_i=0
            for i in range(0,nzz_r_i):
                x1=zb1_r_i[i]
                x2=zb2_r_i[i]
                fx1, dfx1 = ftd_cyl_2_r_i(square,x1,fx1,dfx1,di)
                if (fx1 > 0. and dfx1 < 0.):
                    index=(di.n_kj*3)*i_r_i+3*ii+0
                    d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_i,x1,x2,di.xacc,square,di)
                    i_r_i=i_r_i+1
                elif (fx1 < 0. and dfx1 > 0.):
                    index=(di.n_kj*3)*i_r_i+3*ii+0
                    d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_i,x1,x2,di.xacc,square,di)
                    i_r_i=i_r_i+1
            index=3*ii+0
            d_o.i_bound[index]=i_r_i
#----------------------------------------------------------------
#     Calcolo delle radici immaginarie Caso ref<0 ! kz0 immaginary, kz1 real
        elif (ref < 0.):
            i_i_r=0
            for i in range(0,nzz_i_r):
                x1=zb1_i_r[i]
                x2=zb2_i_r[i]        
                fx1, dfx1 = ftd_cyl_2_i_r(square,x1,fx1,dfx1,di)
                if (fx1 > 0. and dfx1 < 0.):
                    index=(di.n_kj*3)*i_i_r+3*ii+1
                    d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_i_r,x1,x2,di.xacc,square,di)
                    i_i_r=i_i_r+1
                elif (fx1 < 0. and dfx1 > 0.):
                    index=(di.n_kj*3)*i_i_r+3*ii+1
                    d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_i_r,x1,x2,di.xacc,square,di)
                    i_i_r=i_i_r+1
            index=3*ii+1
            d_o.i_bound[index]=i_i_r
#---------------------------------------------------------------
        i_r_r=0 
#-----------------------------------------------------------*
#     classificazione delle discontinuita' per radici reali
        jj=0
        for i in range(1,di.ni+1):
            if (((2*i-1)*pi_2/(di.z0+2*di.A0e/(3*di.ud0))) > kz0min):
                jj=jj+1  
                discontinuita_1[jj-1]=(2*i-1)*pi_2/(di.z0+2*di.A0e/(3*di.ud0))
        iden=0
        i=1
        base=0. 
        while (sqrt(fabs(base)) < discontinuita_1[jj-1]):
            discontinuita_2[i-1]=(2*i-1)*pi_2/(di.z1+2*di.A1e/(3*di.ud1))
            base=((ref*di.n0+pow(discontinuita_2[i-1],2.)*d1_1*di.n0/di.n1)/d0_0)
            if (base > 0.):
                if (sqrt(base) > kz0min):
                    iden=iden+1
                    discontinuita_2[iden-1]=sqrt(base)
            i=i+1
        it=jj+iden	  
#      if (idisc_t == 1)
#      {
#        free(discontinuita_t);    
#	  }
#	  discontinuita_t=calloc(it,sizeof(double));
#      idisc_t=1;
      
        if (idisc_t == 0):
            discontinuita_t=np.zeros(it,dtype=np.float)
            idisc_t=1
        elif(idisc_t == 1):
            del discontinuita_t
            discontinuita_t=np.zeros(it,dtype=np.float)
      
        for i in range(0,jj):
            discontinuita_t[i]=discontinuita_1[i]
	  
        for i in range(0,iden):
            discontinuita_t[jj+i]=discontinuita_2[i]      	  
            
        #selectionSort(it, discontinuita_t);
        discontinuita_t = np.sort(discontinuita_t)

#-----------------------------------------------------------------	  
#     Calcolo delle radici reali kz0 e kz1
        
        nb=6000 
        x1=kz0min+di.precisione
        x2=discontinuita_t[0]-di.precisione	 
        nb = zbrak_cyl_2(ft_cyl_2_r_r,x1,x2,di.nh,xb1,xb2,nb,di,square) 
        if (nb > 0):
            for i in range(0,nb):
                x1=xb1[i]
                x2=xb2[i]
                fx1, dfx1 = ftd_cyl_2_r_r(square,x1,fx1,dfx1,di)
                if (fx1 > 0. and dfx1 < 0.):
                    index=(di.n_kj*3)*i_r_r+3*ii+2
                    d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di)               
                    i_r_r=i_r_r+1 			  
                elif (fx1 < 0. and dfx1 > 0.):
                    index=(di.n_kj*3)*i_r_r+3*ii+2
                    d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di)
                    i_r_r=i_r_r+1			  
			  
        i=1 
        while (i_r_r < di.ni and i <= it-1):
            x1=discontinuita_t[i-1]+di.precisione  
            x2=discontinuita_t[i]-di.precisione;
            fx1, dfx1 = ftd_cyl_2_r_r(square,x1,fx1,dfx1,di)
            fx2, dfx2 = ftd_cyl_2_r_r(square,x2,fx2,dfx2,di)
            if(fx1*fx2 < 0.):
                index=(di.n_kj*3)*i_r_r+3*ii+2
                d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di)               
                i_r_r=i_r_r+1				
            elif (fx1*fx2 > 0.):
                nb=6000
                nb = zbrak_cyl_2(ft_cyl_2_r_r,x1,x2,di.nh,xb1,xb2,nb,di,square)
                if (nb > 0):
                    for ibis in range(0,nb):
                        x1=xb1[ibis]
                        x2=xb2[ibis]
                        fx1, dfx1 = ftd_cyl_2_r_r(square,x1,fx1,dfx1,di)
                        if (fx1 > 0. and dfx1 < 0.):
                            index=(di.n_kj*3)*i_r_r+3*ii+2
                            d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di)
                            i_r_r=i_r_r+1
                        elif(fx1 < 0. and dfx1 > 0.):
                            index=(di.n_kj*3)*i_r_r+3*ii+2
                            d_o.kappa_z0[index]=rtsafe_cyl_2(ftd_cyl_2_r_r,x1,x2,di.xacc,square,di)
                            i_r_r=i_r_r+1
            i=i+1
        index=3*ii+2
        d_o.i_bound[index]=i_r_r
    if (idisc_t == 0):
        del discontinuita_t 	   
# fine funzione f_roots
        
        
def f_plot(di, d_o, dp, i_r):
    #int   i,ii,j,jj,index_a,index_b;
    #double t,kj,kz0;
# questa funzione costruisce il vettore del profilo temporale uttilizzando le radici calcolate dalla funzione f_roots

    for i in range(0,di.n_tpsf):
        t=di.tmin+i*di.dt
        for j in range(0,di.n_kj):
            kj=d_o.kappa_j[j]
            index_b=j*3+0
            for ii in range(0,d_o.i_bound[index_b]):
                index_a=(di.n_kj*3)*ii+3*j+0
                kz0=d_o.kappa_z0[index_a]
                dp.r_tpsf[i][i_r]=dp.r_tpsf[i][i_r]+r_cyl_2_r_i(t,kz0,kj,di)
                #dp.t_tpsf[i][i_r]=dp.t_tpsf[i][i_r]+t_cyl_2_r_i(t,kz0,kj,di) #Commented to speed up execution, since not needed in reflectance configuration program
            index_b=j*3+1
            for ii in range(0,d_o.i_bound[index_b]):
                index_a=(di.n_kj*3)*ii+3*j+1
                kz0=d_o.kappa_z0[index_a]
                dp.r_tpsf[i][i_r]=dp.r_tpsf[i][i_r]+r_cyl_2_i_r(t,kz0,kj,di)
                #dp.t_tpsf[i][i_r]=dp.t_tpsf[i][i_r]+t_cyl_2_i_r(t,kz0,kj,di) #Commented to speed up execution, since not needed in reflectance configuration program
            index_b=j*3+2
            for ii in range(0,d_o.i_bound[index_b]):
                index_a=(di.n_kj*3)*ii+3*j+2
                kz0=d_o.kappa_z0[index_a]
                dp.r_tpsf[i][i_r]=dp.r_tpsf[i][i_r]+r_cyl_2_r_r(t,kz0,kj,di)
                #dp.t_tpsf[i][i_r]=dp.t_tpsf[i][i_r]+t_cyl_2_r_r(t,kz0,kj,di) #Commented to speed up execution, since not needed in reflectance configuration program

# fine funzione f_plot
#-----------------------------------------------------------*
                
def f_plot_Raman(di, did1, did2, dp, dpd1, dpd2, dpR, i_r):
    for i in range(0,di.n_tpsf):
        dpR.t1avg[i][i_r] = (1. - dpd1.r_tpsf[i][i_r] / dp.r_tpsf[i][i_r]) / (di.v0*(did1.ua0 - di.ua0))
        dpR.t2avg[i][i_r] = (1. - dpd2.r_tpsf[i][i_r] / dp.r_tpsf[i][i_r]) / (di.v1*(did2.ua1 - di.ua1))
        dpR.r_tpsf_e[i][i_r] = dp.r_tpsf[i][i_r] * (di.usR0*di.v0*dpR.t1avg[i][i_r] + di.usR1*di.v1*dpR.t2avg[i][i_r])

def selectionSort(n, niz):
    """
    

    Parameters
    ----------
    n : INT
        LENGTH OF ARRAY.
    niz : ARRAY OF DOUBLE
        ARRAY TO BE SORTED IN ASCENDING ORDER.

    Returns
    -------
    None.

    """
    for i in range(0,n-1):
        indmin = i
        for j in range(i+1,n):
            if (niz[j] < niz[indmin]):
                indmin = j
        if (indmin != i):
            temp = niz[i]
            niz[i] = niz[indmin]
            niz[indmin] = temp