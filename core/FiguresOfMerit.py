# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:33:23 2020

@author: ŠUŠNJAR
"""

import tesi_stefan.core.SignalsAndSystems as sis
import numpy as np
import math

PEAK_THRESHOLD = 0.5

def height_of_peak(s_the, s_exp):
    """
    

    Parameters
    ----------
    s_the : SPECTRUM
        THEORETICAL SPECTRUM.
    s_exp : SPECTRUM
        EXPERIMENTALLY RETRIEVED (RECONSTRUCTED) SPECTRUM.

    Returns
    -------
    rel_err_height : FLOAT
        (EXPERIMENTAL PEAK HEIGHT / THEORETICAL PEAK HEIGHT) - 1. (RELATIVE ERROR ON Y-AXIS.)
    rel_err_position : FLOAT
        (EXPERIMENTAL PEAK POSITION / THEORETICAL PEAK POSITION) - 1. (RELATIVE ERROR ON X-AXIS.)

    """
    length = s_the.n_points
    indmax_the = 0
    indmax_exp = 0
    for i in range(1,length):
        if (s_the.s[indmax_the]<s_the.s[i]):
            indmax_the = i
        if (s_exp.s[indmax_exp]<s_exp.s[i]):
            indmax_exp = i
    height_the = s_the.s[indmax_the]
    height_exp = s_exp.s[indmax_exp]
    rel_err_height = height_exp/height_the-1
    rel_err_position = s_exp.l[indmax_exp]/s_exp.l[indmax_the]-1
    return rel_err_height, rel_err_position

def contrast_to_noise(s_the, s_exp):
    length = s_the.n_points
    indmax_the = 0
    indmax_exp = 0
    noise = sis.Spectrum('Noise = Experimental - Theoretical')
    noise.l = np.zeros(s_the.n_points,dtype=np.float)
    noise.s = np.zeros(s_the.n_points,dtype=np.float)
    var = 0
    for i in range(1,length):
        if (s_the.s[indmax_the]<s_the.s[i]):
            indmax_the = i
        if (s_exp.s[indmax_exp]<s_exp.s[i]):
            indmax_exp = i
        noise.l[i] = s_the.l[i]
        noise.s[i] = s_exp.s[i]-s_the.s[i]
    i = indmax_the
    while(s_the.s[i]>=s_the.s[indmax_the]*PEAK_THRESHOLD):
        var = var + noise.s[i]**2
        i = i-1
    left = i+1
    i = indmax_the+1
    while(s_the.s[i]>=s_the.s[indmax_the]*PEAK_THRESHOLD):
        var = var + noise.s[i]**2
        i = i+1
    right = i-1
    var = var/(right-left+1)
    #noise.plot("Raman shift [1/cm]","noise")
    return s_exp.s[indmax_exp]/math.sqrt(var)

def suppression(s1_the, s1_exp, s2_the, s2_exp):
    length = s1_the.n_points
    dl = s1_the.l[1]-s1_the.l[0]
    indmax1_the = 0
    indmax1_exp = 0
    indmax2_the = 0
    indmax2_exp = 0
    for i in range(1,length):
        if (s1_the.s[indmax1_the]<s1_the.s[i]):
            indmax1_the = i
        if (s1_exp.s[indmax1_exp]<s1_exp.s[i]):
            indmax1_exp = i
        if (s2_the.s[indmax2_the]<s2_the.s[i]):
            indmax2_the = i
        if (s2_exp.s[indmax2_exp]<s2_exp.s[i]):
            indmax2_exp = i
    A1 = 0
    A2 = 0
    i1 = indmax1_the
    while(s1_the.s[i1]>=s1_the.s[indmax1_the]*PEAK_THRESHOLD):
        A1 = A1+abs(s2_exp.s[i1])*dl
        i1 = i1-1
    i1 = indmax1_the+1
    while(s1_the.s[i1]>=s1_the.s[indmax1_the]*PEAK_THRESHOLD):
        A1 = A1+abs(s2_exp.s[i1])*dl
        i1 = i1+1
    i2 = indmax2_the
    while(s2_the.s[i2]>=s2_the.s[indmax2_the]*PEAK_THRESHOLD):
        A2 = A2+abs(s2_exp.s[i2])*dl
        i2 = i2-1
    i2 = indmax2_the+1
    while(s2_the.s[i2]>=s2_the.s[indmax2_the]*PEAK_THRESHOLD):
        A2 = A2+abs(s2_exp.s[i2])*dl
        i2 = i2+1
    return A2/A1
