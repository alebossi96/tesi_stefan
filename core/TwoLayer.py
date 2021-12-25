# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:55:37 2020

@author: ŠUŠNJAR
"""
import core.Functions as fn
from core.DataStructures import* 
import copy

di = DataInput("file_input.txt")
d_o = DataOutputRaw(di)

did1 = copy.deepcopy(di)
did2 = copy.deepcopy(di)
did1.ua0 = 1.0001*di.ua0 #relative increment of absorption coefficient is 10^(-4) for numerical derivative calculation
did2.ua1 = 1.0001*di.ua1 #relative increment of absorption coefficient is 10^(-4) for numerical derivative calculation
d_o_d1 = copy.deepcopy(d_o)
d_o_d2 = copy.deepcopy(d_o)

dp = DataPlot()
dpd1 = DataPlot()
dpd2 = DataPlot()
dpR = DataPlotRaman()

fn.f_roots(di,d_o)
fn.f_roots(did1,d_o_d1)
fn.f_roots(did2,d_o_d2)

for i_r in range(0,di.n_rec):
    di.ro = di.rec[i_r]
    did1.ro = did1.rec[i_r]
    did2.ro = did2.rec[i_r]
    fn.f_plot(di,d_o,dp,i_r)
    fn.f_plot(did1,d_o_d1,dpd1,i_r)
    fn.f_plot(did2,d_o_d2,dpd2,i_r)
    fn.f_plot_Raman(di,did1,did2,dp,dpd1,dpd2,dpR,i_r)

if (di.i_type==3):
    dpR.print_to_file(di)
else:
    dp.print_to_file(di)


