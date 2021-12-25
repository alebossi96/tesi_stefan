# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:03:36 2020

@author: ŠUŠNJAR
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import core.ResolutionToPoints as rtp
class Spectrum(object):
    #lambda - l
    #value of spectrum - s
    #number of points - n_points
    def __init__(self,info, instrumentData):
        """
        

        Parameters
        ----------
        info : STRING
            INFORMATION ABOUT THE SPECTRUM.

        Returns
        -------
        None.

        """
        self.info = info
        self.instrumentData = instrumentData
    
    def read_from_file(self,filename,separator=' '):
        """
        EACH LINE OF FILE SHOULD BE IN FORMAT '<lambda><separator><value>'. 
        WAVELENGTH AND VALUE ARE SEPARATED WITH 1 SEPARATOR CHARACTER (ASSUMED SPACE BY DEFAULT). 
        ACCORDING TO THE NUMBER OF LINES IN THE FILE, TWO ARRAYS FOR WAVELENGTHS AND VALUES ARE DEFINED AND FILLED WITH THE READ VALUES.

        Parameters
        ----------
        filename : STRING
            NAME OF FILE TO BE READ.
        separator : CHAR
            SEPARATING CHARACTER BETWEEN WAVELENGTH VALUE AND SPECTRUM VALUE.

        Returns
        -------
        None.

        """
        self.instrumentData.n_points  = 0
        f1 = open(filename,"r")
        while(True):
            if (f1.readline()==""):
                break
            else:
                self.instrumentData.n_points +=1
        f1.close()
        f2 = open(filename,"r")
        self.l = np.zeros(self.instrumentData.n_points ,dtype=np.float)
        self.s = np.zeros(self.instrumentData.n_points ,dtype=np.float)
        for i in range(0, self.instrumentData.n_points ):
            string = f2.readline()
            ind = string.find(separator)
            self.l[i] = (float)(string[:ind])
            self.s[i] = (float)(string[ind+1:-1])
        f2.close()
    
    def set_gaussian(self, mean, amplitude):
        """
        

        Parameters
        ----------
        mean : FLOAT
            CENTRAL WAVELENGTH.
        sigma : FLOAT
            STANDARD DEVIATION.
        number_of_points : INT
            NUMBER OF POINTS IN SPECTRUM - THE HIGHER THE BETTER RESOLUTION.
        lambda_min : FLOAT
            SMALLEST WAVELENGTH IN STORED SPECTRUM.
        lambda_max : FLOAT
            LARGEST WAVELENGTH IN STORED SPECTRUM.
        amplitude : FLOAT
            IN GAUSSIAN DISTRIBUTION S(lambda) = A*exp(-(lambda-mean)^2/(2*sigma^2), FACTOR A.

        Returns
        -------
        None.

        """

        if self.instrumentData.n_points  == 1:
            self.l = np.zeros(1,dtype=np.float)
            self.s = np.zeros(1,dtype=np.float)
            self.l[0] = mean
            self.s[0] = amplitude
        else:
            max_wn = self.instrumentData.max_wavenumber
            min_wn = self.instrumentData.min_wavenumber
            dl = (max_wn-min_wn)/(self.instrumentData.n_points -1)
            self.l = np.zeros(self.instrumentData.n_points ,dtype=np.float)
            self.s = np.zeros(self.instrumentData.n_points ,dtype=np.float)
            for i in range(0,self.instrumentData.n_points ):
                wn = min_wn+i*dl
                self.l[i] = wn
                self.s[i] = amplitude*math.exp(-pow(wn-mean,2)/(2*pow(self.instrumentData.sigma,2)))
    
    def copy_from(self,spec,limit=-1):
        """
        

        Parameters
        ----------
        spec : SPECTRUM
            THE SPECTRUM WHICH WAVELENGTHS AND VALUES ARE TO BE COPIED.
        limit : INT
            MAXIMUM INDEX UP TO WHICH THE COPYING IS DONE. IF limit==-1, COPYING IS DONE UNTIL THE END OF spec ARRAYS. Default is -1.

        Returns
        -------
        None.

        """
        
        if (limit>=0 and limit<spec.n_points):
            self.instrumentData.n_points  = limit+1
            self.l = np.zeros(self.instrumentData.n_points ,dtype=np.float)
            self.s = np.zeros(self.instrumentData.n_points ,dtype=np.float)
            for i in range(0,self.instrumentData.n_points ):
                self.l[i] = spec.l[i]
                self.s[i] = spec.s[i]
        else:
            self.instrumentData.n_points  = spec.n_points
            self.l = np.zeros(self.instrumentData.n_points ,dtype=np.float)
            self.s = np.zeros(self.instrumentData.n_points ,dtype=np.float)
            for i in range(0,self.instrumentData.n_points ):
                self.l[i] = spec.l[i]
                self.s[i] = spec.s[i]
            if (limit != -1):
                print("Unexpected value for limit. Index is outside of bounds. The whole spectrum is copied.")
        
        
    def fill_with_gaussian(self, mean, sigma, amplitude):
        """
        

        Parameters
        ----------
        mean : FLOAT
            CENTRAL WAVELENGTH.
        sigma : FLOAT
            STANDARD DEVIATION.
        amplitude : FLOAT
            IN GAUSSIAN DISTRIBUTION S(lambda) = A*exp(-(lambda-mean)^2/(2*sigma^2), FACTOR A.

        Returns
        -------
        None.

        """
        for i in range(0,self.instrumentData.n_points ):
                lam = self.l[i]
                self.s[i] = amplitude*math.exp(-pow(lam-mean,2)/(2*pow(sigma,2)))
        
    def add_spectrum(self,S2,info=""):
        """
        

        Parameters
        ----------
        S2 : SPECTRUM
            SPECTRUM TO BE ADDED.
        info : STRING, optional
            IF CHANGE IN info IS NEEDED, PLEASE GIVE THE VALUE TO THIS ARGUMENT DIFFERENT FROM "". The default is "".

        Returns
        -------
        None.

        """
        
        S1new, S2new = merge_spectra(self, S2)
        self.instrumentData.n_points  = S1new.n_points
        if (info != ""):
            self.info = info
        for i in range(0,self.instrumentData.n_points ):
            S1new.s[i] = S1new.s[i]+S2new.s[i]
        self.l = S1new.l
        self.s = S1new.s
        
    def scale(self,factor):
        """
        

        Parameters
        ----------
        factor : FLOAT
            SCALING FACTOR BY WHICH ALL THE VALUES OF SPECTRUM ARE MULTIPLIED.

        Returns
        -------
        None.

        """
        
        for i in range(0,self.instrumentData.n_points ):
            self.s[i] *= factor
        
    
    def plot(self,xlabel,ylabel):
        plt.title(self.info) 
        plt.xlabel(xlabel) 
        plt.ylabel(ylabel)
        plt.plot(self.l,self.s,scaley='log')
        plt.show()

    def first_derivative(self):
        fd = Spectrum(self.info+' - first derivative')
        fd.n_points = self.instrumentData.n_points 
        fd.l = np.zeros(self.instrumentData.n_points ,dtype = np.float)
        fd.s = np.zeros(self.instrumentData.n_points ,dtype = np.float)
        for i in range(0, self.instrumentData.n_points -1):
            fd.l[i] = self.l[i]
            fd.s[i] = (self.s[i+1]-self.s[i])/(self.l[i+1]-self.l[i])
        fd.l[self.instrumentData.n_points -1] = self.l[self.instrumentData.n_points -1]
        fd.s[self.instrumentData.n_points -1] = 0
        return fd
        

class Signal(object):
    #s - value array (1D if n_points==0, 2D if n_points>0)
    #t - time 1D array
    #l - wavelength 1D array (only defined if n_points>0)
    #n_points - number of wavelength array points (keep 0 if signal is not decomposed by wavelengths)
    
    def __init__(self,length,n_points=0):
        """
        

        Parameters
        ----------
        n_points : INT
            NUMBER OF DIFFERENT WAVELENGTH COMPONENTS OF SIGNAL. USE 0 (DEFAULT) IF THE SIGNAL IS NOT DECOMPOSED IN WAVELENGTHS.
        length : INT
            NUMBER OF DISCRETE TIME INSTANTS REPRESENTING THE SIGNAL IN DIGITAL DOMAIN.

        Returns
        -------
        None.

        """
        self.instrumentData.n_points  = n_points
        if (n_points > 0):
            self.s = np.zeros((n_points,length),dtype=np.float)
            self.l = np.zeros(n_points,dtype=np.float)
        else:
            self.s = np.zeros(length,dtype=np.float)
        self.t = np.zeros(length,dtype=np.float)
    
    def read_from_file(self,filename,separator=' '):
        """
        EACH LINE OF FILE SHOULD BE IN FORMAT '<time><separator><value>'. 
        TIME AND VALUE ARE SEPARATED WITH 1 SEPARATOR CHARACTER (ASSUMED SPACE BY DEFAULT). 
        ACCORDING TO THE NUMBER OF LINES IN THE FILE, TWO ARRAYS FOR TIMES AND VALUES ARE DEFINED AND FILLED WITH THE READ VALUES.

        Parameters
        ----------
        filename : STRING
            NAME OF FILE TO BE READ.
        separator : CHAR
            SEPARATING CHARACTER BETWEEN TIME VALUE AND INSTRUMENT RESPONSE FUNCTION VALUE.

        Returns
        -------
        None.

        """
        length = 0
        f1 = open(filename,"r")
        while(True):
            if (f1.readline()==""):
                break
            else:
                length+=1
        f1.close()
        f2 = open(filename,"r")
        self.instrumentData.n_points  = 0 #Signal is not decomposed by wavelengths
        self.t = np.zeros(length,dtype=np.float)
        self.s = np.zeros(length,dtype=np.float)
        for i in range(0, length):
            string = f2.readline()
            ind = string.find(separator)
            self.t[i] = (float)(string[:ind])
            self.s[i] = (float)(string[ind+1:-1])
        f2.close()
    
    def read_from_file_decomposed(self,filename,separator=','):
        """
        THE FIRST LINE OF FILE CONTAINS ARBITRARY (LABELLING) FIELD, THEN <separator><l1><separator><l2>...<separator><ln>.
        FROM THE SECOND, EACH LINE OF FILE SHOULD BE IN FORMAT '<time><separator><v1><separator><v2><separator>...<vn>'. 
        <v1>...<vn> are values at <time> at wavelengths <l1>...<ln>
        ACCORDING TO THE NUMBER OF LINES IN THE FILE AND THE NUMBER OF VALUES IN THE FIRST ROW, THREE ARRAYS FOR TIMES, WAVELENGHTS AND VALUES ARE DEFINED AND FILLED WITH THE READ VALUES.

        Parameters
        ----------
        filename : STRING
            NAME OF FILE TO BE READ.
        separator : CHAR
            SEPARATING CHARACTER BETWEEN VALUES IN THE SAME ROW.

        Returns
        -------
        None.

        """
        
        f1 = open(filename,"r")
        string = f1.readline()
        ind = string.find(separator)
        n_points = 0
        while(ind>0):
            n_points+=1
            string = string[ind+1:]
            ind = string.find(separator)
        length = 0
        while(True):
            if (f1.readline()==""):
                break
            else:
                length+=1
        f1.close()
        f2 = open(filename,"r")
        self.instrumentData.n_points  = n_points
        self.t = np.zeros(length,dtype=np.float)
        self.l = np.zeros(n_points,dtype=np.float)
        self.s = np.zeros((n_points,length),dtype=np.float)
        string = f2.readline()
        ind = string.find(separator)
        string = string[ind+1:]
        for j in range(0, n_points-1):
            ind = string.find(separator)
            self.l[j] = (float)(string[:ind])
            string = string[ind+1:]
        self.l[n_points-1] = (float)(string[:-1])
        for i in range(0, length):
            string = f2.readline()
            ind = string.find(separator)
            self.t[i] = (float)(string[:ind])
            string = string[ind+1:]
            for j in range(0, n_points-1):
                ind = string.find(separator)
                self.s[j][i] = (float)(string[:ind])
                string = string[ind+1:]
            self.s[n_points-1][i] = (float)(string[:-1])
        f2.close()
    
    def plot_at_wavelength(self,wavelength,title,xlabel,ylabel):
        if (self.instrumentData.n_points  > 0):
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            ind = find_nearest_wavelength(self.l,wavelength)
            plt.plot(self.t,self.s[ind])
            plt.show()
        else:
            print("Error. Signal is not prepared in a proper format. Try plotting using the function: plot.")
    
    def plot(self,title,xlabel,ylabel):
        if (self.instrumentData.n_points  == 0):
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.plot(self.t,self.s)
            plt.show()
        else:
            print("Error. Signal is not prepared in a proper format. Try plotting at a specific wavelength, since this signal is decomposed by wavelengths. Use the function: plot_at_wavelength.")

class System(object):
    def __init__(self):
        self.irf = Signal(1)
    
    def read_from_file(self,filename,separator=' '):
        self.irf.read_from_file(filename,separator)
    
    def set_gaussian_irf(self, mean, sigma, number_of_points, t_min, t_max, amplitude, normalize=True):
        """
        

        Parameters
        ----------
        mean : float
            CENTRAL TIME INSTANT (MAXIMUM OF IRF).
        sigma : float
            STANDARD DEVIATION.
        number_of_points : int
            NUMBER OF POINTS (TIME INSTANTS) IN IRF.
        t_min : float
            SMALLEST TIME INSTANT STORED IN IRF.
        t_max : float
            LARGEST TIME INSTANT STORED IN IRF.
        amplitude : float
            IN GAUSSIAN DISTRIBUTION IRF(t) = A*exp(-(t-mean)^2/(2*sigma^2), FACTOR A.
        normalize : BOOLEAN, optional
            IF True, VALUES ARE NORMALIZED TO GIVE TOTAL SUM (INTEGRAL) EQUAL TO amplitude; IF False, VALUES ARE DETERMINED ACCORDING TO THE FORMULA GIVEN ABOVE, WITHOUT ANY NORMALIZATION. The default is True.

        Returns
        -------
        None.

        """
        suma = 0
        if (number_of_points <= 0 or not(number_of_points==(int)(number_of_points))):
            print("Error. Positive integer expected for number of points.")
        elif (number_of_points == 1):
            self.irf.n_points = 0
            self.irf.t = np.zeros(number_of_points,dtype=np.float)
            self.irf.s = np.zeros(number_of_points,dtype=np.float)
            self.irf.t[0] = mean
            self.irf.s[0] = amplitude
        else:
            self.irf.n_points = 0
            dt = (t_max-t_min)/(number_of_points-1)
            self.irf.t = np.zeros(number_of_points,dtype=np.float)
            self.irf.s = np.zeros(number_of_points,dtype=np.float)
            for i in range(0,number_of_points):
                t = t_min+i*dt
                self.irf.t[i] = t
                self.irf.s[i] = amplitude*math.exp(-pow(t-mean,2)/(2*pow(sigma,2)))
                suma += self.irf.s[i]
        
        if (normalize and number_of_points>1):
            k = amplitude/suma #scaling factor
            for i in range(0, number_of_points):
                self.irf.s[i] *= k
           
            
    
    def filt(self, signal):
        """
        

        Parameters
        ----------
        s : SIGNAL
            INPUT SIGNAL. SAMPLING PERIODS OF s AND IMPULSE RESPONSE FUNCTION OF THE SYSTEM irf SHOULD BE EQUAL.

        Returns
        -------
        out : SIGNAL
            OUTPUT SIGNAL, OBTAINED AFTER PASSING THE INPUT SIGNAL THROUGH THE FILTER DEFINED BY THE IMPULSE RESPONSE FUNCTION. IMPULSE RESPONSE FUNCTION OF THIS FILTER IS DEFINED BY THE INSTRUMENT RESPONSE FUNCTION OF THE SYSTEM.

        """
        s = signal
        ind0 = find_nearest_wavelength(self.irf.t, s.t[0]) #The same method as for wavelengths, but here we are looking for times
        indmax = len(self.irf.t)-1-ind0
        indmin = -ind0
        #Here we assume that time instants in signal s and in IRF of the filter are equidistant (same sampling period)
        #Otherwise, this method will not work properly
        len_s = len(s.t)
        len_irf = len(self.irf.t)
        len_out = len_s+len_irf-1
        out = Signal(len_out,s.n_points)
        dt = s.t[1]-s.t[0]
        # if (s.n_points>0):
        #     for l in range(0,s.n_points):
        #         out.l[l] = s.l[l]
        #         for n in range(indmin, len_s+indmax):
        #             out.s[l][n-indmin] = 0
        #             out.t[n-indmin] = s.t[0] + n*dt
        #             for k in range(0,len_s):
        #                 if (n-k>=indmin and n-k<=indmax):
        #                     out.s[l][n-indmin] += s.s[l][k]*self.irf.s[n-k]
        # else:
        #     out.s[n-indmin] = 0
        #     out.t[n-indmin] = s.t[0] + n*dt
        #     for k in range(0,len_s):
        #         if (n-k>=indmin and n-k<=indmax):
        #             out.s[n-indmin] += s[k]*self.irf.s[n-k]
        if (s.n_points>0):
            for l in range(0,s.n_points):
                out.l[l] = s.l[l]
                for n in range(0, len_out):
                    out.s[l][n] = 0
                    out.t[n] = s.t[0] + self.irf.t[0] + n*dt
                    for k in range(0,len_s):
                        if (n-k>=0 and n-k<=len_irf-1):
                            out.s[l][n] += s.s[l][k]*self.irf.s[n-k]
        else:
            for n in range(0,len_out):
                out.s[n] = 0
                out.t[n] = s.t[0] + self.irf.t[0] + n*dt
                for k in range(0,len_s):
                    if (n-k>=0 and n-k<=len_irf-1):
                        out.s[n] += s[k]*self.irf.s[n-k]
        
        return out
    
    def gated_signal(self, signal, n_gates, tmin, tmax):
        """
        

        Parameters
        ----------
        s : SIGNAL
            INPUT SIGNAL WHICH IS TO BE GATED IN TIME.
        n_gates : INT
            NUMBER OF TIME GATES.
        tmin : FLOAT
            START OF THE FIRST GATE.
        tmax : FLOAT
            END OF THE LAST GATE.

        Returns
        -------
        out : SIGNAL
            GATED SIGNAL, LENGTH IN TIME IS EQUAL TO THE NUMBER OF GATES n_gates.

        """
        s = signal
        if (tmin<=s.t[0] or tmax>=s.t[-1] or tmin>=tmax):
            print("Error. Start and end time of gates should be inside the interval on which the signal s is defined. Input signal s is returned.")
            return s
        dt = s.t[1]-s.t[0]
        dt_gate = (tmax-tmin)/n_gates #Duration of single time gate
        out = Signal(n_gates, s.n_points)
        if (s.n_points>0):
            for l in range(0, s.n_points):
                i = 0
                m = 0
                out.l[l] = s.l[l]
                while(s.t[i]<tmin):
                    i+=1
                out.t[0] = tmin
                out.s[l][0] = (s.t[i]-tmin)*s.s[l][i-1]
                while(True):
                    while(s.t[i]<(m+1)*dt_gate+tmin):
                        out.s[l][m] += s.s[l][i]*dt
                        i += 1
                    delta = s.s[l][i-1]*(s.t[i]-(m+1)*dt_gate-tmin) #This has to be written in the next gate and subtracted from the current gate
                    out.s[l][m] -= delta
                    m+=1
                    if (m==n_gates):
                        break
                    out.s[l][m] = delta
                    out.t[m] = tmin + m*dt_gate
        else:
            i = 0
            m = 0
            while(s.t[i]<tmin):
                i+=1
            out.t[0] = tmin
            out.s[0] = (s.t[i]-tmin)*s.s[i-1]
            while(True):
                while(s.t[i]<(m+1)*dt_gate+tmin):
                    out.s[m] += s.s[i]*dt
                    i += 1
                delta = s.s[i-1]*(s.t[i]-(m+1)*dt_gate-tmin) #This has to be written in the next gate and subtracted from the current gate
                out.s[m] -= delta
                m+=1
                if (m==n_gates):
                    break
                out.s[m] = delta
                out.t[m] = tmin + m*dt_gate
                
        return out
    
    def add_Poisson_noise(self, counts,signal):
        """
        

        Parameters
        ----------
        counts : INT
            TOTAL NUMBER OF PHOTON COUNTS.
        signal : SIGNAL
            SIGNAL TO WHICH THE (POISSON) NOISE IS TO BE ADDED.

        Returns
        -------
        out : SIGNAL
            INPUT SIGNAL + (POISSON) NOISE.
        gain : FLOAT
            CALCULATED GAIN NEEDED TO OBTAIN TOTAL NUMBER OF COLLECTED PHOTONS EQUAL TO counts.

        """
        total = 0
        length = len(signal.t)
        out = Signal(length,signal.n_points)
        gain = 0
        if (signal.n_points == 0):
            for i in range(0, length):
                total = total + signal.s[i] #Should we add here *deltaT?
            gain = counts/total
            for k in range(0, length):
                out.s[k] = np.random.poisson(abs(signal.s[k]*gain))
                out.t[k] = signal.t[k]
        else:
            for i in range(0, length):
                for j in range(0, signal.n_points):
                    total = total + signal.s[j][i] #Should we add here *deltaT?
            gain = counts/total
            print(gain)
            print(total)
            print(counts)
            for l in range(0, signal.n_points):
                out.l[l] = signal.l[l]
                for k in range(0, length):
                    out.s[l][k] = np.random.poisson(abs(signal.s[l][k]*gain)) #abs added to take care of eventual negative spectral measurements
                    out.t[k] = signal.t[k]
        return out, gain
        
                
    
    
def find_nearest_wavelength(lambda_array, wavelength):
    """
    

    Parameters
    ----------
    lambda_array : ARRAY OF FLOAT
        ARRAY OF WAVELENGTHS.
    wavelength : FLOAT
        DESIRED WAVELENGTH.

    Returns
    -------
    ind : INT
        INDEX OF ELEMENT OF lambda_array WHICH IS EQUAL (OR THE NEAREST ONE) TO THE SPECIFIED wavelength. IN CASE OF ERRORS, -1 IS RETURNED.

    """
    
    n = len(lambda_array)
    i=0
    while(lambda_array[i]<wavelength):
        i+=1
        if (i>=n):
            return -1
    lambdaUpper = lambda_array[i]
    if (lambdaUpper==wavelength):
        return i
    else:
        if (i==0): 
            return i
        elif (wavelength-lambda_array[i-1]<lambdaUpper-wavelength):
            return i-1
        else:
            return i
    return -1


def merge_spectra(S1, S2):
    """
    

    Parameters
    ----------
    S1 : SPECTRUM
        SPECTRUM 1.
    S2 : SPECTRUM
        SPECTRUM 2.

    Returns
    -------
    S1new : SPECTRUM
        NEW SPECTRUM 1, EXTENDED ON WAVELENGTHS CONTAINED IN SPECTRUM 2 THAT ARE MISSING IN SPECTRUM 1.
    S2new : SPECTRUM
        NEW SPECTRUM 2, EXTENDED ON WAVELENGTHS CONTAINED IN SPECTRUM 1 THAT ARE MISSING IN SPECTRUM 2.

    """
    length = S1.n_points+S2.n_points
    S1new = Spectrum(S1.info)
    S2new = Spectrum(S2.info)
    S1new.n_points = length
    S2new.n_points = length
    S1new.l = np.zeros(length,dtype=np.float)
    S1new.s = np.zeros(length,dtype=np.float)
    S2new.l = np.zeros(length,dtype=np.float)
    S2new.s = np.zeros(length,dtype=np.float)
    i = 0
    j = 0
    k = 0
    real_length = length
    while (k<real_length):
        if  (S1.l[i]<S2.l[j]):
            S1new.l[k] = S1.l[i]
            S2new.l[k] = S1.l[i]
            S1new.s[k] = S1.s[i]
            if (j==0):
                S2new.s[k] = 0
            else:
                S2new.s[k] = S2.s[j-1]*(S2.l[j]-S1.l[i])/(S2.l[j]-S2.l[j-1])+S2.s[j]*(S1.l[i]-S2.l[j-1])/(S2.l[j]-S2.l[j-1])
            i += 1
        elif (S2.l[j]<S1.l[i]):
            S1new.l[k] = S2.l[j]
            S2new.l[k] = S2.l[j]
            S2new.s[k] = S2.s[j]
            if (i==0):
                S1new.s[k] = 0
            else:
                S1new.s[k] = S1.s[i-1]*(S1.l[i]-S2.l[j])/(S1.l[i]-S1.l[i-1])+S1.s[i]*(S2.l[j]-S1.l[i-1])/(S1.l[i]-S1.l[i-1])
            j += 1
        else: #They are equal
            S1new.l[k] = S1.l[i]
            S2new.l[k] = S2.l[j]
            S1new.s[k] = S1.s[i]
            S2new.s[k] = S2.s[j]
            i += 1
            j += 1
            real_length -= 1
        k += 1
        if (i==S1.n_points):
            while(k<real_length):
                S1new.l[k] = S2.l[j]
                S2new.l[k] = S2.l[j]
                S1new.s[k] = 0
                S2new.s[k] = S2.s[j]
                j += 1
                k += 1
            break
        if (j==S2.n_points):
            while(k<real_length):
                S1new.l[k] = S1.l[i]
                S2new.l[k] = S1.l[i]
                S1new.s[k] = S1.s[i]
                S2new.s[k] = 0
                i += 1
                k += 1
            break
    if (real_length<length):
        S1new_cut = Spectrum(S1new.info)
        S2new_cut = Spectrum(S2new.info)
        S1new_cut.copy_from(S1new,real_length-1)
        S2new_cut.copy_from(S2new,real_length-1)
        return S1new_cut, S2new_cut
    return S1new, S2new

def total_number_of_collected_photons(sig, continuous = False):
    """
    

    Parameters
    ----------
    sig : SIGNAL
        SIGNAL THAT REPRESENTS THE NUMBER OF PHOTONS COLLECTED DURING TIME (IN CASE OF A REAL MEASUREMENT), OR SIGNAL THAT REPRESENTS SIMULATED, THEORETICALLY CALCULATED SIGNAL AT CERTAIN TIME INSTANTS (IN LATTER CASE, NUMERICAL INTEGRATION IS PERFORMED WHEN CALCULATING TOTAL NUMBER OF COLLECTED PHOTONS).
    continuous : BOOL, optional
        True IF SIGNAL IS THEORETICALLY OBTAINED (SIMULATED), BY CALCULATING IN DIFFERENT TIME INSTANTS, UNIT IS 1/TIME, OR (1/TIME)/AREA. False IF SIGNAL IS EXPERIMENTALLY OBTAINED, BY MEASURING THE NUMBER OF PHOTONS REACHING THE DETECTOR (IT IS BY NATURE DIMENSIONLESS SIGNAL, OR EXPRESSED AS 1/AREA). The default is False.

    Returns
    -------
    suma : INT
        NUMBER OF COLLECTED PHOTONS (EITHER BY REAL MEASUREMENT, OR OBTAINED NUMERICALLY AFTER INTEGRATION OVER TIME).

    """
    
    suma = 0
    length = len(sig.t)
    for i in range(0,sig.n_points):
        for j in range(0,length-1):
            if (continuous):
                suma += sig.s[i][j] * (sig.t[j+1]-sig.t[j]) #Multiplying by deltaT - integration
            else:
                suma += sig.s[i][j] #Only discrete sum (for real measurements)
        if (not continuous):
            suma += sig.s[i][length-1]
    return suma
    
