# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 14:53:55 2020

@author: ŠUŠNJAR
"""
import tesi_stefan.core.Functions as fn
from tesi_stefan.core.DataStructures import*
import copy
import numpy as np
from tesi_stefan.core.SignalsAndSystems import*
import os

class ForwardModel(object):
    
    def __init__(self, geometry, time_step = 15.03, n_tpsf =256,
                 input_file = os.path.dirname(__file__) + '/file_input.txt', relative_increment_absorption = 0.0001, update=False):
        """
        

        Parameters
        ----------
        input_file : STRING, optional
            NAME OF THE INPUT FILE THAT CONTAINS ALL THE INFORMATION ABOUT THE MEDIA. The default is 'file_input.txt'.
        relative_increment_absorption : FLOAT, optional
            USED FOR COMPUTATION OF NUMERICAL DERIVATIVES REQUIRED FOR CALCULATION OF AVERAGE TIMES t1 AND t2 THAT LIGHT 'SPENDS' IN TWO LAYERS. The default is 0.0001.
        update : BOOLEAN, optional
            True IF UPDATE IS NEEDED (VALUES USED FOR CALCULATION DIFFER FROM THOSE IN THE INPUT FILE). PROGRAMMER IS RESPONSIBLE TO RUN THE METHOD ForwardModel.calculate() IN THIS CASE (IF update == True), AFTER OVERRIDING DEFAULT PARAMETERS. IN CASE update == False, NOTHING IS NEEDED, EVERYTHING IS DONE AUTOMATICALLY. The default is False.

        Returns
        -------
        None.

        """
        #TODO parte 1,2 critiche
        
        self.di = DataInput(input_file, geometry, time_step, n_tpsf)
        self.relative_increment_absorption = relative_increment_absorption
        
        if (not update):
            self.d_o = DataOutputRaw(self.di)
            self.did1 = copy.deepcopy(self.di)
            self.did2 = copy.deepcopy(self.di)
            self.did1.ua0 = (1+relative_increment_absorption)*self.di.ua0 #relative increment of absorption coefficient is 10^(-4) for numerical derivative calculation
            self.did2.ua1 = (1+relative_increment_absorption)*self.di.ua1 #relative increment of absorption coefficient is 10^(-4) for numerical derivative calculation
            self.d_o_d1 = copy.deepcopy(self.d_o)
            self.d_o_d2 = copy.deepcopy(self.d_o)
            
            self.dp = DataPlot()
            self.dpd1 = DataPlot()
            self.dpd2 = DataPlot()
            self.dpR = DataPlotRaman()
            #TODO critical
            #TODO capire se queste funzioni esistono già fatte
            fn.f_roots(self.di,self.d_o)
            fn.f_roots(self.did1,self.d_o_d1)
            fn.f_roots(self.did2,self.d_o_d2)
            
            for i_r in range(0,self.di.n_rec):
                self.di.ro = self.di.rec[i_r]
                self.did1.ro = self.did1.rec[i_r]
                self.did2.ro = self.did2.rec[i_r]
                fn.f_plot(self.di,self.d_o,self.dp,i_r)
                fn.f_plot(self.did1,self.d_o_d1,self.dpd1,i_r)
                fn.f_plot(self.did2,self.d_o_d2,self.dpd2,i_r)
                fn.f_plot_Raman(self.di,self.did1,self.did2,self.dp,self.dpd1,self.dpd2,self.dpR,i_r)

            #Calculation of sensitivity matrix W(x)
            self.W1 = np.zeros(self.di.n_tpsf,dtype=np.float)
            self.W2 = np.zeros(self.di.n_tpsf,dtype=np.float)
            for i in range(0,self.di.n_tpsf):
                self.W1[i] = self.dp.r_tpsf[i][0]*self.di.v0*self.dpR.t1avg[i][0]
                self.W2[i] = self.dp.r_tpsf[i][0]*self.di.v1*self.dpR.t2avg[i][0]
    def calculate(self):
        self.d_o = DataOutputRaw(self.di)
        self.did1 = copy.deepcopy(self.di)
        self.did2 = copy.deepcopy(self.di)
        self.did1.ua0 = (1+self.relative_increment_absorption)*self.di.ua0 #relative increment of absorption coefficient is 10^(-4) for numerical derivative calculation
        self.did2.ua1 = (1+self.relative_increment_absorption)*self.di.ua1 #relative increment of absorption coefficient is 10^(-4) for numerical derivative calculation
        self.d_o_d1 = copy.deepcopy(self.d_o)
        self.d_o_d2 = copy.deepcopy(self.d_o)
        
        self.dp = DataPlot()
        self.dpd1 = DataPlot()
        self.dpd2 = DataPlot()
        self.dpR = DataPlotRaman()
        
        fn.f_roots(self.di,self.d_o)
        fn.f_roots(self.did1,self.d_o_d1)
        fn.f_roots(self.did2,self.d_o_d2)
        
        for i_r in range(0,self.di.n_rec):
            self.di.ro = self.di.rec[i_r]
            self.did1.ro = self.did1.rec[i_r]
            self.did2.ro = self.did2.rec[i_r]
            fn.f_plot(self.di,self.d_o,self.dp,i_r)
            fn.f_plot(self.did1,self.d_o_d1,self.dpd1,i_r)
            fn.f_plot(self.did2,self.d_o_d2,self.dpd2,i_r)
            fn.f_plot_Raman(self.di,self.did1,self.did2,self.dp,self.dpd1,self.dpd2,self.dpR,i_r)
        
        
        #Calculation of sensitivity matrix W(x)
        self.W1 = np.zeros(self.di.n_tpsf,dtype=np.float)
        self.W2 = np.zeros(self.di.n_tpsf,dtype=np.float)
        for i in range(0,self.di.n_tpsf):
            self.W1[i] = self.dp.r_tpsf[i][0]*self.di.v0*self.dpR.t1avg[i][0]
            self.W2[i] = self.dp.r_tpsf[i][0]*self.di.v1*self.dpR.t2avg[i][0]
    
    def sensitivity_matrix(self):
        """
        IF RAMAN REFLECTANCE R(lambda,x) = W(x) * [S1(lambda); S2(lambda)], WHERE S1(lambda) AND S2(lambda) ARE SPECTRA OF THE TWO LAYERS, AND W(x) IS THE SENSITIVITY MATRIX, THEN W1(x) AND W2(x) ARE THE 2 COLUMNS OF THIS MATRIX.

        Returns
        -------
        W1 : ARRAY OF FLOAT
            SENSITIVITY MATRIX (ACTUALLY COLUMN VECTOR) FOR THE FIRST LAYER. 
        W2 : ARRAY OF FLOAT
            SENSITIVITY MATRIX (ACTUALLY COLUMN VECTOR) FOR THE SECOND LAYER.

        """
        
        W1 = np.zeros(self.di.n_tpsf,dtype=np.float)
        W2 = np.zeros(self.di.n_tpsf,dtype=np.float)
        for i in range(0,self.di.n_tpsf):
            W1[i] = self.dp.r_tpsf[i][0]*self.di.v0*self.dpR.t1avg[i][0]
            W2[i] = self.dp.r_tpsf[i][0]*self.di.v1*self.dpR.t2avg[i][0]
        return W1, W2
    #Method sensitivity_matrix is not recommended for calling. 
    #Values of W1 and W2 are already calculated during initialization and stored as self.W1 and self.W2 of the object of class ForwardModel
    
    def output(self,Spectrum_Top_layer, Spectrum_Bottom_layer):
        """
        

        Parameters
        ----------
        S1 : SPECTRUM
            RAMAN SPECTRUM OF THE FIRST LAYER.
        S2 : SPECTRUM
            RAMAN SPECTRUM OF THE SECOND LAYER (SAME WAVELENGTH ARRAYS ARE EXPECTED FOR S2 AND S1).

        Returns
        -------
        sig : SIGNAL
            THEORETICALLY CALCULATED RAMAN REFLECTANCE SIGNAL.

        """
        S1 = Spectrum_Top_layer
        S2 = Spectrum_Bottom_layer
        sig = Signal(self.di.n_tpsf,S1.instrumentData.n_points)
        for i in range(0,S1.instrumentData.n_points):
            sig.l[i] = S1.l[i]
            for j in range(0,self.di.n_tpsf):
                sig.t[j] = self.di.tmin+j*self.di.dt
                sig.s[i][j] = self.W1[j]*S1.s[i] + self.W2[j]*S2.s[i]
        
        return sig
    
    def gated_sensitivity_matrix(self, n_gates, tmin, tmax):
        """
        

        Parameters
        ----------
        n_gates : INT
            NUMBER OF TIME GATES.
        tmin : FLOAT
            START OF THE FIRST GATE.
        tmax : FLOAT
            END OF THE LAST GATE.

        Returns
        -------
        Wg1 : ARRAY OF FLOAT
            GATED SENSITIVITY MATRIX (ACTUALLY COLUMN VECTOR) FOR THE FIRST LAYER. LENGTH IS EQUAL TO THE NUMBER OF GATES n_gates.
        Wg2 : ARRAY OF FLOAT
            GATED SENSITIVITY MATRIX (ACTUALLY COLUMN VECTOR) FOR THE SECOND LAYER. LENGTH IS EQUAL TO THE NUMBER OF GATES n_gates.

        """
        dt = self.di.dt
        if (tmin<=self.di.tmin or tmax>=self.di.tmin+(self.di.n_tpsf-1)*dt or tmin>=tmax):
            print(tmax,self.di.tmin+(self.di.n_tpsf-1)*dt)
            print(tmax)
            print("Error. Start and end time of gates should be inside the interval on which the model output is defined. Unchanged self.W1 and self.W2 are returned from this object of class ForwardModel.")
            return self.W1, self.W2
        
        dt_gate = (tmax-tmin)/n_gates #Duration of single time gate
        Wg1 = np.zeros(n_gates,dtype=np.float)
        Wg2 = np.zeros(n_gates,dtype=np.float)
        i = 0
        m = 0
        while(self.di.tmin+i*dt<tmin):
            i+=1
        Wg1[0] = (self.di.tmin+i*dt-tmin)*self.W1[i-1]
        Wg2[0] = (self.di.tmin+i*dt-tmin)*self.W2[i-1]
        while(True):
            while(self.di.tmin+i*dt<(m+1)*dt_gate+tmin):
                Wg1[m] += self.W1[i]*dt
                Wg2[m] += self.W2[i]*dt
                i += 1
            delta1 = self.W1[i-1]*(self.di.tmin+i*dt-(m+1)*dt_gate-tmin) #This has to be written in the next gate and subtracted from the current gate
            Wg1[m] -= delta1
            delta2 = self.W2[i-1]*(self.di.tmin+i*dt-(m+1)*dt_gate-tmin) #This has to be written in the next gate and subtracted from the current gate
            Wg2[m] -= delta2
            m+=1
            if (m==n_gates):
                break
            Wg1[m] = delta1
            Wg2[m] = delta2
        
        return Wg1, Wg2
        
        
        
