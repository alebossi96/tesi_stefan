# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:11:02 2020

@author: ŠUŠNJAR
"""

errorFormat = "Unexpected format of input file!\n"
import core.Functions as fn
import math
import numpy as np

N_REC_MAX = 7
N_PUNTI_TPSF_MAX = 400


class DataInput(object):
    #z0 = 5 # thickness of the first layer (mm)
    #z1 = 17 # thickness of the second layer (mm)
    ##rec[N_REC_MAX] = 10;  # set of source-receiver distances
    #rec = [10, 10, 10, 10, 10, 10, 10]
    #ro = 10     # point where the TPSF is evaluated (distance from the laser beam mm)
    #R = 50      # lateral dimention of the cylinder (Radius of the cylinder mm)
    #ua0 = 0.011    # absorption of the first layer (mm^-1) (only background uab0)
    #ua1 = 0.003    # absorption of the second layer (mm^-1) (only background uab1)
    #ud0 = 1.65    # reduced scattering coefficient of the first layer (mm^-1) 
    #ud1 = 1.65    # reduced scattering coefficient of the second layer (mm^-1)
    #usR0 = 0.002    # Raman scattering coefficient of the first layer (mm^-1)
    #usR1 = 0.002    # Raman scattering coefficient of the second layer (mm^-1)
    #n0 = 1.4     # refractive index of the first layer 
    #n1 = 1.4     # refractive index of the second layer 
    #ne = 1.0     # refractive index of the external medium 
    #c = 0.299792458      # speed of light in air (mm/ps)
    #v0 = c/n0     # speed of light in the first layer (mm/ps)
    #v1 = c/n1     # speed of light in the second layer (mm/ps) 
    #n0e = n0/ne   # relative refractive index first layer - external medium
    #n1e = n1/ne   # relative refractive index second layer - external medium
    #an01 = 1   # relative refractive index first-second layer
    #A0e = 1    # A factor between first layer and external
    #A1e = 1	 # A factor between second layer and external
    #dt = 20     # time step of the TPSF  (ps)
    #tmin = 100	 # Starting time of the TPSF  (ps)
    #dr = 1
    #rmin = 1
    #xacc =  1e-10 # accuracy of the roots calculated with the routine rtsafe
    #precisione = 1e-10; # accuracy
    #pi = 1     # 2*asin(1) wrong value here obviously
    #i_type = 3   # reflectance (1) or transmittance (2) or Raman refl. (3)
    #n_rec = 1  # number of receiver
    #n_tpsf = 300  # number of steps in the TPSF
    #n_cw = 0
    #ni = 35     # number of roots kn0 and kn1 
    #n_kj = 35 # number of term kj 
    #nh = 5000     # partition number of the interval for real roots 
    #nhi = 500    # partition number of the interval for immaginary roots
    #outputFile = "out.csv" #output file name
    def __init__(self, filename):
        f = open(filename, "r")
        key = "nome_file"
        string = ""
        c_ext =".csv";
        self.rec = np.zeros(N_REC_MAX,dtype=np.float)
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1): 
            print(errorFormat)
        else:
            self.outputFile = string+c_ext
   
        key = "indice_analisi"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.i_type = int(string)
        
        key = "spessore_0"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.z0 = float(string)
            
        key = "spessore_1"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.z1 = float(string)

        key = "numero_ricevitori"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.n_rec = int(string)
   
        key = "distanza_1"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rec[0] = float(string)
   
        key = "distanza_2"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rec[1] = float(string)
      
        key = "distanza_3"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rec[2] = float(string)
        
        key = "distanza_4"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rec[3] = float(string)
        
        key = "distanza_5"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rec[4] = float(string)
            
        key = "distanza_6"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rec[5] = float(string)
            
        key = "distanza_7"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rec[6] = float(string)
            
        key = "dimensione_laterale"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.R = float(string)
   
        key = "assorbimento_0"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.ua0 = float(string)
 
        key = "assorbimento_1"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.ua1 = float(string)
   
        key = "scattering_0"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.ud0 = float(string)   
      
        key = "scattering_1"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.ud1 = float(string)

        key = "Raman_scattering_0"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.usR0 = float(string)
      
        key = "Raman_scattering_1"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.usR1 = float(string)
   
        key = "indice_r_0"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.n0 = float(string)    
      
        key = "indice_r_1"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.n1 = float(string)
   
        key = "indice_r_e"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.ne = float(string)    

        key = "velocita_vuoto"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.c = float(string)
        
        key = "time_step"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.dt = float(string)    
   
        key = "time_min"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.tmin = float(string)
   
        key = "distance_step"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.dr = float(string)   
   
        key = "distance_min"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.rmin = float(string)
   
        key = "accuratezza_roots"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.xacc = float(string)   

        key = "accuratezza_range"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.precisione = float(string)   
   
        key = "punti_tpsf"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.n_tpsf = int(string)
   
        key = "numero_roots_kn0"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.ni = int(string)
   
        key = "numero_kj"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.n_kj = int(string)
   
        key = "partizione_reale"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.nh = int(string)   
   
        key = "partizione_immaginaria"
        status, string = fn.read_data_line_from_file(f, key)
        if (status == -1):
            print(errorFormat)
        else:
            self.nhi = int(string)    

        self.n0e=self.n0/self.ne;
        self.n1e=self.n1/self.ne;
        self.A0e=fn.dAA_cyl(self.n0e);
        self.A1e=fn.dAA_cyl(self.n1e);
        self.pi=2*math.asin(1.);
        self.v0=self.c/self.n0;
        self.v1=self.c/self.n1;
        
        f.close()
import os
cwd = os.getcwd()
class DataOutputRaw(object):
    # Structure of data output for the program cyl_lay2_turbo C version
	# up dated on june 2003
	# Version modified for Python
    #kappa_z0 = []     # vettore dove vengono immagazzinate le radici kz0
    #kappa_j = []		# vettore dove vengono immagazzianate le radici kj
    #i_bound = []      # vettore che contiene il numero di radici kz0 per ogni kj
    
    def __init__(self, di):
        self.kappa_j = np.zeros(di.n_kj,dtype=np.float)
        self.i_bound = np.zeros(di.n_kj*3,dtype=np.int)
        self.kappa_z0 = np.zeros(di.n_kj*3*di.ni,dtype=np.float)
        import os
        f=open(os.path.dirname(__file__)+"/j0.txt","r")
        for i in range(0,di.n_kj):
            self.kappa_j[i] = ((float)(f.readline()))/di.R
        f.close()
        
class DataPlot(object):
	# Structure of data output for the program cyl_lay2_turbo C version
	# up dated on june 2003
	# Version modified for Python

    #r_tpsf = np.zeros((N_PUNTI_TPSF_MAX,N_REC_MAX),dtype=np.float)  # vettore dove vengono immagazzinate la riflettanza risolta in tempo
    #t_tpsf = np.zeros((N_PUNTI_TPSF_MAX,N_REC_MAX),dtype=np.float)  # vettore dove vengono immagazzinate la trasmittanza risolta in tempo
	
    def __init__(self):
        self.r_tpsf = np.zeros((N_PUNTI_TPSF_MAX,N_REC_MAX),dtype=np.float)
        self.t_tpsf = np.zeros((N_PUNTI_TPSF_MAX,N_REC_MAX),dtype=np.float)
    
    def print_to_file(self,di):
        f = open(di.outputFile,"w")
        f.write("Characteristics of the medium\n")
        if (di.i_type==1):
            f.write("Reflectance data\n")
        elif (di.i_type==2):
            f.write("Transmittance data\n")
        else:
            return
        f.write("%s,%f\n" % ("Lateral dimension [mm]=", di.R))
        f.write("%s,%f\n" % ("Thickness first layer [mm]=", di.z0))
        f.write("%s,%f\n" % ("Thickness second layer [mm]=", di.z1))
        f.write("%s,%f\n" % ("Absorption first layer [mm^-1]=", di.ua0))
        f.write("%s,%f\n" % ("Absorption second layer [mm^-1]=", di.ua1))
        f.write("%s,%f\n" % ("Reduced scattering first layer [mm^-1]=", di.ud0))
        f.write("%s,%f\n" % ("Reduced scattering second layer [mm^-1]=", di.ud1))
        f.write("%s,%f\n" % ("Raman scattering first layer [mm^-1]=", di.usR0))
        f.write("%s,%f\n" % ("Raman scattering second layer [mm^-1]=", di.usR1))
        f.write("%s,%f\n" % ("Refraction index first layer=", di.n0))
        f.write("%s,%f\n" % ("Refraction index second layer=", di.n1))
        f.write("%s,%f\n" % ("Refraction index external medium=", di.ne))
        f.write("%s,%i\n" % ("Nnumber of receivers=", di.n_rec))
        f.write("Distances considered [mm]=,")
        for i in range(0, di.n_rec-1):
            f.write("%f," % di.rec[i])
        f.write("%f\n" % di.rec[di.n_rec-1])
        
        for i in range(0,di.n_tpsf):
            t = di.tmin+i*di.dt
            for i_r in range(0,di.n_rec-1):
                if (di.i_type == 1):
                    f.write("%f, %e" % (t, self.r_tpsf[i][i_r]))
                elif (di.i_type == 2):
                    f.write("%f, %e" % (t, self.t_tpsf[i][i_r]))
            i_r = di.n_rec-1
            if (di.i_type == 1):
                f.write("%f, %e\n" % (t, self.r_tpsf[i][i_r]))
            elif (di.i_type == 2):
                f.write("%f, %e\n" % (t, self.r_tpsf[i][i_r]))
        
        f.close()
        
class DataPlotRaman(object):
	# Class for Raman reflectance data output
	# 2020 Python version

    #r_tpsf_e = [[0]*N_REC_MAX]*N_PUNTI_TPSF_MAX  # array where the fluence (at Raman emission wavelength) vs time is stored for the reflectance configuration
    #t1avg = [[0]*N_REC_MAX]*N_PUNTI_TPSF_MAX  # array where the average time spent by detected light in the first layer is stored
    #t2avg = [[0]*N_REC_MAX]*N_PUNTI_TPSF_MAX  # array where the average time spent by detected light in the second layer is stored
	
    def __init__(self):
        self.r_tpsf_e = np.zeros((N_PUNTI_TPSF_MAX,N_REC_MAX),dtype=np.float)
        self.t1avg = np.zeros((N_PUNTI_TPSF_MAX,N_REC_MAX),dtype=np.float)
        self.t2avg = np.zeros((N_PUNTI_TPSF_MAX,N_REC_MAX),dtype=np.float)
    
    def print_to_file(self,di):
        f = open(di.outputFile,"w")
        f.write("Characteristics of the medium\n")
        if (di.i_type==3):
            f.write("Reflectance Raman data\n")
        else:
            return
        f.write("%s,%f\n" % ("Lateral dimension [mm]=", di.R))
        f.write("%s,%f\n" % ("Thickness first layer [mm]=", di.z0))
        f.write("%s,%f\n" % ("Thickness second layer [mm]=", di.z1))
        f.write("%s,%f\n" % ("Absorption first layer [mm^-1]=", di.ua0))
        f.write("%s,%f\n" % ("Absorption second layer [mm^-1]=", di.ua1))
        f.write("%s,%f\n" % ("Reduced scattering first layer [mm^-1]=", di.ud0))
        f.write("%s,%f\n" % ("Reduced scattering second layer [mm^-1]=", di.ud1))
        f.write("%s,%f\n" % ("Raman scattering first layer [mm^-1]=", di.usR0))
        f.write("%s,%f\n" % ("Raman scattering second layer [mm^-1]=", di.usR1))
        f.write("%s,%f\n" % ("Refraction index first layer=", di.n0))
        f.write("%s,%f\n" % ("Refraction index second layer=", di.n1))
        f.write("%s,%f\n" % ("Refraction index external medium=", di.ne))
        f.write("%s,%i\n" % ("Nnumber of receivers=", di.n_rec))
        f.write("Distances considered [mm]=,")
        for i in range(0, di.n_rec-1):
            f.write("%f," % di.rec[i])
        f.write("%f\n" % di.rec[di.n_rec-1])
        f.write("t[ps],Raman reflectance[mm^{-2}ps^{-1}],<t1>[ps],<t2>[ps] \n")
        
        for i in range(0,di.n_tpsf):
            t = di.tmin+i*di.dt
            for i_r in range(0,di.n_rec-1):
                f.write("%f, %e, %e, %e" % (t, self.r_tpsf_e[i][i_r], self.t1avg[i][i_r], self.t2avg[i][i_r]))
            i_r = di.n_rec-1
            f.write("%f, %e, %e, %e\n" % (t, self.r_tpsf_e[i][i_r], self.t1avg[i][i_r], self.t2avg[i][i_r]))
                   
        f.close()
        
        
        
        
        
        
        
