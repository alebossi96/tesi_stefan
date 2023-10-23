import numpy as np
import math 
import scipy
import scipy.special
import core.ForwardModel as fwd
import core.Geometry as geo
import core_f as core
import matplotlib.pyplot as plt
def geometry(mu_s_top, mu_s_bottom):
    return geo.Geometry(thickness_top = 10, thickness_bottom = 10,
                                        mu_a_top = 0.0075, mu_a_bottom = 0.002,
                                        mu_s_top = mu_s_top, mu_s_bottom = mu_s_bottom,
                                        mu_R_top = 0.00011, mu_R_bottom = 0.011,
                                        r0=4.5, r1 = 5.5, n_step_r = 2)
def convolution_gaus_decay(t, x):
    #TODO essendoci il delay di 60ps da problemi?
    #no aggiunta costanti
    #sig = sig ^2
    mu = x[0]
    sig = x[1]
    tau = x[2] #tempo rilassamento
    res = np.exp(-t/tau)*(1-scipy.special.erf((sig+2*tau*(mu-t))/(2*tau*np.sqrt(sig))))
    return res
def res_irf(t, geo, top, Wg1, x):
    C_M = core.conv_matrix(convolution_gaus_decay(t,x), len(t))
    rec = np.matmul(C_M, Wg1)
    return top/np.max(top)- rec/np.max(rec)
def get_irf(t, geo, top, Wg):
    #(Wg1, Wg2) = get_matrix_el(geo, t)
    residual = lambda x: res_irf(t, geo, top, Wg, x)
    x0 = np.array([t[np.argmax(top)],8000, 157])
    #plt.plot(Wg2/np.max(Wg2))
    return scipy.optimize.least_squares(fun = residual, x0 = x0)
def get_param_from_irf(t, data):
    residual = lambda x: convolution_gaus_decay(t, x)
    x0 = np.array([t[np.argmax(data)],8000, 157])
    return scipy.optimize.least_squares(fun = residual, x0 = x0)
    
def reconstruction(x, t, geom, M_T, M_B):
    #print(x)
    mu = x[0]
    sig = x[1]
    tau = x[2] #tempo rilassamento
    len_wl = int((len(x)-3)/2)
    mu_R_T = x[3:(len_wl+3)]
    mu_R_B = x[(len_wl+3):(2*len_wl+3)] 
    G_T = np.matmul(M_T, convolution_gaus_decay(t, [mu, sig, tau]))
    G_B = np.matmul(M_B, convolution_gaus_decay(t, [mu, sig, tau]))
    
    #return  mult_mu_R_G(mu_R_T, G_T)+mult_mu_R_G(mu_R_B, G_B)
    res = np.matmul(mu_R_T[:,np.newaxis], G_T[np.newaxis,:])+np.matmul(mu_R_B[:,np.newaxis], G_B[np.newaxis,:])
    return res/np.sum(res)  
    
def get_param_spectra(time, data, geom):
    forward = fwd.ForwardModel(geom, time_step = time[1]-time[0], n_tpsf = len(time))
    (W_1, W_2) = forward.sensitivity_matrix()
    #TODO W_1, W_2 vanno calcolati solo 1 volta!!!!
    M_T = core.conv_matrix(W_1, len(time))
    M_B = core.conv_matrix(W_2, len(time))
    norm_dif = lambda x: np.linalg.norm(reconstruction(x, time, geom, M_T, M_B)-data/np.sum(data))
    num_wl = (len(data),)
    x0 = np.zeros((3+len(data)+len(data),))
    x0[0] = time[np.argmax(np.sum(data,axis = 0))]
    x0[1] = 8000
    x0[2] = 157
    (W_1, _) = get_matrix_el(geom, time)
    sum_time = np.sum(W_1)
    
    init_mu_R = np.sum(data, axis = 1)/sum_time#TODO fare meglio init mu_R
    #print(init_mu_R)
    for i in range(len(init_mu_R)):
        x0[i+3] = init_mu_R[i]
        x0[i+3+len(data)] = init_mu_R[i]
    #x0 = np.array([time[np.argmax(inp_data)],8000, 0.4, 1.2, np.ones(num_wl), np.ones(num_wl)])#è un np array!
    lower_bound = [0]*(3+2*len(data))
    upper_bound = [np.inf]*(3+2*len(data))
    bounds = (lower_bound, upper_bound) 
    return scipy.optimize.least_squares(fun = norm_dif, x0 = x0)#,  max_nfev = 40, diff_step = 0.01,  bounds = bounds,  verbose = 2
def conv_gaus_decay_green(t, x, layer):
    #TODO essendoci il delay di 60ps da problemi?
    #no aggiunta costanti
    #sig = sig ^2
    #print(x)
    mu = x[0]
    sig = x[1]
    tau = 157 #tempo rilassamento
    if x[2]> 0:
        mu_s_top =  x[2] 
    else:
        mu_s_top = 10
    
    mu_s_top = 0.4
    if x[3]> 0:
        mu_s_bottom =  x[3] 
    else:
        mu_s_bottom = 0.1
    geom = geometry(mu_s_top, mu_s_bottom)
    (W_1, W_2) = get_matrix_el(geom, t)
    if layer == "B":
        M_W = core.conv_matrix(W_2, len(t))
    elif layer == "T":
        M_W = core.conv_matrix(W_1, len(t))
    res = np.matmul(M_W, convolution_gaus_decay(t, [mu, sig]))
    return res/np.max(res)
def ideal_conv_irf_green_fun(t, geometry, data):
    (W_1, _) = get_matrix_el(geometry, t)
    M_W = core.conv_matrix(W_1, len(t))
    fun = lambda t, x: np.matmul(M_W,convolution_gaus_decay(t,x))    
    return fit_param(fun, data, t)
    
def get_irf_green_fun_mu_s(t, data_top, data_bottom ):  
    return fit_param_mu_s(conv_gaus_decay_green, data_top, data_bottom, t)
    
    
def fit_param_mu_s(fun, data_top, data_bottom, time):
    fun_t = lambda x: fun(time, x, "T")
    fun_b = lambda x: fun(time, x, "B")
    residual = lambda x: 2*np.abs(fun_t(x) - data_top/np.max(data_top))+ np.abs(fun_b(x) - data_bottom/np.max(data_bottom))
    x0 = np.array([time[np.argmax(data_top)]-100,1e4, 0.4, 1.2])
    return scipy.optimize.least_squares(fun = residual, x0 = x0, max_nfev = 30, jac = '3-point', verbose = 2)
    
#TODO fare recons qui
def mult_mu_R_G(mu_R, G):
    mat = np.zeros((len(mu_R),len(G)))
    for i in range(len(mu_R)):
        mat[i,:] = mu_R[i]
    return np.matmul(mat, G)

def fit_param(fun, inp_data, time):
    residual = lambda x: fun(time, x) - inp_data/np.max(inp_data)
    x0 = np.array([time[np.argmax(inp_data)],8000])
    return scipy.optimize.least_squares(fun = residual, x0 = x0)
    
def get_matrix_el(geometry, time):
    W_t_1 = np.zeros((len(time),))
    W_t_2 = np.zeros((len(time),))
    dt = time[1]-time[0]
    for r in geometry:
        forward = fwd.ForwardModel(geometry, time_step = dt, n_tpsf = len(time))
        Wg1, Wg2 = forward.sensitivity_matrix()
        W_t_1 += r*Wg1*2*math.pi*geometry.dr
        W_t_2 += r*Wg2*2*math.pi*geometry.dr
    return (W_t_1,  W_t_2)    
