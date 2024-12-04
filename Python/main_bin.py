"""  Run DAC simulation using various noise shaping techniques

@author: Bikash Adhikari 
@date: 23.02.2024
@license: BSD 3-Clause
"""

# %%
import numpy as np
from scipy  import signal
import scipy
import csv
import matplotlib.pyplot as plt
import statistics
import math 
import tqdm
from static_dac_model import quantise_signal, generate_code, generate_dac_output
from  quantiser_configurations import quantiser_configurations, get_measured_levels

from MHOQ import MHOQ
from MHOQ_BIN import MHOQ_BIN
from nsdcal import nsdcal, noise_shaping
# from dsm import dsm
from filter_parameters import filter_params

# from plot_variance import plot_variance, plot_enob
from process_sim_output import process_sim_output
from balreal import balreal
from welch_psd import welch_psd
from simulation_plots import bar_plot 
from figures_of_merit import TS_SINAD


# Generate test signal
def test_signal(SCALE, MAXAMP, FREQ, Rng,  OFFSET, t):
    Xcs = (SCALE/100)*MAXAMP*np.cos(2*np.pi*FREQ*t) + OFFSET  +Rng/2
    return Xcs 

# %% Choose MHOQ implemetatio method
class MHOQ_IMPLEMENTATION:
    BINARY = 1
    SCALED = 2
# MHOQ_METHOD = MHOQ_IMPLEMENTATION.BINARY
MHOQ_METHOD = MHOQ_IMPLEMENTATION.SCALED

# %% Quantiser configurations 
Qconfig =1 
Nb, Mq, Vmin, Vmax, Rng, Qstep, YQ, Qtype = quantiser_configurations(Qconfig, MHOQ_METHOD)

# %% Sampling frequency and rate
Fs = 1e6
Ts = 1/Fs

# %% Generate time vector
Xcs_SCALE = 100
Xcs_FREQ = 999
match 2:
    case 1:  # specify duration as number of samples and find number of periods
        Nts = 1e5  # no. of time samples
        Np = np.ceil(Xcs_FREQ*Ts*Nts).astype(int) # no. of periods for carrier

    case 2:  # specify duration as number of periods of carrier
        Np = 10 # no. of periods for carrier
        
Npt = 1  # no. of carrier periods to use to account for transients
Np = Np + Npt

t_end = Np/Xcs_FREQ  # time vector duration
t = np.arange(0, t_end, Ts)  # time vector

# %% Generate carrier/test signal
SIGNAL_MAXAMP = Rng/2
SIGNAL_OFFSET = -Qstep/2  # try to center given quantiser type
Xcs = test_signal(Xcs_SCALE, SIGNAL_MAXAMP, Xcs_FREQ, Rng,  SIGNAL_OFFSET, t)

# %% Reconstruction filter
N_lp = 2        # Filter order
FcS = 1e5       # Filter cutoff frequency
b_r, a_r = signal.butter(N_lp, FcS/(Fs/2),"low")

# %% Choose filter type for filter parameters to use in MHOQ impelmentation
class RCF: # 
    Perception = 1          # Psychoacoustically optimal noise-shaping filter form Goodwin et. al 2003
    Butterworth = 2         # Butterworth filter
    Butterworth_H2_Hinif = 3 # ||NTF||_{2} <= 1.5           # Optimal filter obtained minimising the error variance. 
    Butterworth_H2_Hinif_Ext = 4 # ||NTF||_{2} <= 1.5           # Optimal filter obtained minimising the error variance. 
    Butterworth_H2_H2 = 5 # ||R[z]-1||<= 1.5^2           # Optimal filter obtained minimising the error variance. 
    Butterworth_H2_H2 = 6 # ||R[z]-1||<= 2.5           # Optimal filter obtained minimising the error variance. 
    Butterworth_H2_H2_Ext = 7 # ||R[z]-1||<= 1.5^2           # Optimal filter obtained minimising the error variance. 
    Butterworth_H2_H2_Ext = 8 # ||R[z]-1||<= 2.5           # Optimal filter obtained minimising the error variance. 
# % Get  filter parameters - Choose the options from above.
LPF = 3   # Low pass filter type 
b_lpf, a_lpf, Fc = filter_params(LPF, FcS)

# %% State space matrices
# Tf to state space
A, B, C, D  = signal.tf2ss(b_lpf, a_lpf)
A, B, C, D = balreal(A, B, C, D) # balanced realisation
# eigenvalues of the state matrix A
eig_A, eig_V = np.linalg.eig(A)

# %% Quatniser Model
# Quantiser model: 1 - Ideal , 2- Calibrated
QMODEL = 1

# %% Quantiser levels : Ideal and Measured
YQns = YQ.squeeze()
INL  =  get_measured_levels(Qconfig).squeeze()
MLns = YQns + INL*Qstep # Measured quantiser levels without measurement error

# Add measurement error to INL
sd_inl = 0.01*statistics.stdev(INL)
Measr_err = np.round(np.random.uniform(-sd_inl,sd_inl,len(INL)),4)
INL_ME = INL + Measr_err
# MLns_ME = YQns + INL_ME*Qstep   # Measured quantiser levels with measurement error
MLns_ME = MLns

# %% LIN methods on/off
DIR_ON = True
NSD_ON = True
MPC_ON = True

# %% Direct Quantization 
if DIR_ON:
    # C_DIR = direct_quant(Xcs, Qstep, Vmin, Vmax)
    C_DIR = quantise_signal(Xcs, Qstep, YQns, Qtype)
    match QMODEL:
        case 1:
            Xcs_DIR = generate_dac_output(C_DIR, YQns)
        case 2:
            Xcs_DIR = generate_dac_output(C_DIR, MLns_ME)

# %% NSD CAL
if NSD_ON:
    b_nsf = b_lpf - a_lpf
    a_nsf = b_lpf
    C_NSD = noise_shaping(Nb, Xcs, b_nsf, a_nsf, Qstep, YQns, MLns_ME, Vmin, QMODEL)  
    match QMODEL:
        case 1:
            Xcs_NSD = generate_dac_output(C_NSD, YQns)
        case 2:
            Xcs_NSD = generate_dac_output(C_NSD, MLns) 

# # %% MPC : Prediction horizon
match MHOQ_METHOD:
    case 1:
        MHOQ = MHOQ_BIN(Nb, Qstep, QMODEL, A, B, C, D)
    case 2:
        MHOQ = MHOQ(Nb, Qstep, QMODEL, A, B, C, D)

N1 = 1
C_MHOQ1, ObjVal1 = MHOQ.get_codes(N1, Xcs, YQ, MLns_ME)
match QMODEL:
    case 1:
        Xcs_MHOQ1 = generate_dac_output(C_MHOQ1, YQ)
    case 2:
        Xcs_MHOQ1 = generate_dac_output(C_MHOQ1, MLns) 

N2 = 2
C_MHOQ2, ObjVal2  = MHOQ.get_codes(N2, Xcs, YQ, MLns_ME)
match QMODEL:
    case 1:
        Xcs_MHOQ2 = generate_dac_output(C_MHOQ2, YQ)
    case 2:
        Xcs_MHOQ2 = generate_dac_output(C_MHOQ2, MLns) 


# N3 = 3
# C_MHOQ3, ObjVal3  = MHOQ.get_codes(N3, Xcs, YQ, MLns_ME)
# match QMODEL:
#     case 1:
#         Xcs_MHOQ3 = generate_dac_output(C_MHOQ3, YQ)
#     case 2:
#         Xcs_MHOQ3 = generate_dac_output(C_MHOQ3, MLns) 


# N4 = 4
# C_MHOQ4, ObjVal4  = MHOQ.get_codes(N4, Xcs, YQ, MLns_ME)
# match QMODEL:
#     case 1:
#         Xcs_MHOQ4 = generate_dac_output(C_MHOQ4, YQ)
#     case 2:
#         Xcs_MHOQ4 = generate_dac_output(C_MHOQ4, MLns) 


# %% Trim vector lengths  
tm = t[:Xcs_MHOQ2.shape[1]]
len_tm = len(tm)
TRANSOFF = np.floor(Npt*Fs/Xcs_FREQ).astype(int)  # remove transient effects from output

# %% 

SINAD_COMP_SEL = 1
plot_val = False
FXcs_DQ, SINAD_DQ, ENOB_DQ = process_sim_output(tm, Xcs_DIR, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot= plot_val, descr="DQ")
FXcs_NSD, SINAD_NSD, ENOB_NSD= process_sim_output(tm, Xcs_NSD, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="NSD")
FXcs_MHOQ1, SINAD_MHOQ1, ENOB_MHOQ1= process_sim_output(tm, Xcs_MHOQ1, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ1")
FXcs_MHOQ2, SINAD_MHOQ2, ENOB_MHOQ2= process_sim_output(tm, Xcs_MHOQ2, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ2")
# FXcs_MHOQ3, SINAD_MHOQ3, ENOB_MHOQ3= process_sim_output(tm, Xcs_MHOQ3, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ3")
# FXcs_MHOQ4, SINAD_MHOQ4, ENOB_MHOQ4= process_sim_output(tm, Xcs_MHOQ4, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ4")
# %% Simulation plots
# from simulation_plots import plot_variance
# bar_plot(descr= "ENOB", DQ =  ENOB_DQ, NSD = ENOB_NSD, MHOQ1 = ENOB_MHOQ1, MHOQ2 = ENOB_MHOQ2, MHOQ3 = ENOB_MHOQ3, MHOQ4 = ENOB_MHOQ4 )
# bar_plot(descr= "SINAD", DQ =  SINAD_DQ, NSD = SINAD_NSD, MHOQ1 = SINAD_MHOQ1, MHOQ2 = SINAD_MHOQ2, MHOQ3 = SINAD_MHOQ3, MHOQ4 = SINAD_MHOQ4 )


bar_plot(descr= "ENOB", DQ =  ENOB_DQ, NSD = ENOB_NSD, MHOQ1 = ENOB_MHOQ1, MHOQ2 = ENOB_MHOQ2)
bar_plot(descr= "SINAD", DQ =  SINAD_DQ, NSD = SINAD_NSD, MHOQ1 = SINAD_MHOQ1, MHOQ2 = SINAD_MHOQ2)

# %% Error- Unfiltered  
err_DIR = Xcs - Xcs_DIR.squeeze()
err_NSD = Xcs - Xcs_NSD.squeeze()
err_MHOQ1 = Xcs.squeeze()[:len(Xcs_MHOQ1.squeeze())] - Xcs_MHOQ1.squeeze()
err_MHOQ2 = Xcs[:len(Xcs_MHOQ2.squeeze())] - Xcs_MHOQ2.squeeze()

err_DIR = err_DIR[TRANSOFF:-TRANSOFF]
err_NSD = err_NSD[TRANSOFF:-TRANSOFF]
err_MHOQ1 = err_MHOQ1[TRANSOFF:-TRANSOFF]
err_MHOQ2 = err_MHOQ2[TRANSOFF:-TRANSOFF]
# %%


var_DIR = 1/(len(err_DIR))*sum(i**2 for i in err_DIR)
var_NSD = 1/(len(err_NSD))*sum(i**2 for i in err_NSD)
var_MHOQ1 = 1/(len(err_MHOQ1))*sum(i**2 for i in err_MHOQ1)
var_MHOQ2 = 1/(len(err_MHOQ2))*sum(i**2 for i in err_MHOQ2)


# %% 
# sl = 1000
# fig,ax = plt.subplots()
# ax.plot(t[:sl], Xcs_DIR.squeeze()[:sl])
# ax.plot(t[:sl], Xcs_NSD.squeeze()[:sl])


# %%
# brecos = b_r
# arecos = a_r
# F_Xcs_DIR = signal.lfilter(brecos, arecos, Xcs_DIR)
# F_Xcs_NSD = signal.lfilter(brecos, arecos, Xcs_NSD)
# F_Xcs_MHOQ1 = signal.lfilter(brecos, arecos, Xcs_MHOQ1)
# F_Xcs_MHOQ2 = signal.lfilter(brecos, arecos, Xcs_MHOQ2)

# SINAD_DQ = TS_SINAD(F_Xcs_DIR.squeeze()[TRANSOFF:-TRANSOFF], t[TRANSOFF:-TRANSOFF], plot_val, 'DQ')
# SINAD_NSD = TS_SINAD(F_Xcs_NSD.squeeze()[TRANSOFF:-TRANSOFF], t[TRANSOFF:-TRANSOFF], plot_val, 'NSD')
# SINAD_MHOQ1 = TS_SINAD(F_Xcs_MHOQ1.squeeze()[TRANSOFF:-TRANSOFF], t[TRANSOFF:-TRANSOFF-1], plot_val, 'MHOQ1')
# SINAD_MHOQ2 = TS_SINAD(F_Xcs_MHOQ2.squeeze()[TRANSOFF:-TRANSOFF], t[TRANSOFF:-TRANSOFF-2], plot_val, 'MHOQ2')

# bar_plot(descr= "SINAD", DQ =  SINAD_DQ, NSD = SINAD_NSD, MHOQ1 = SINAD_MHOQ1, MHOQ2 = SINAD_MHOQ2)
# %% frequency response 
w_ntf, h_ntf = signal.freqz(a_lpf, b_lpf, fs= Fs)  #w - angular frequency, h= frequency response
h_ntf_db = 20*np.log10(h_ntf)

# %% PSD plots of the given signal 
# f, Pxx = signal.periodogram(err_DIR, Fs)
# Pxx = 20*np.log10(Pxx)

# # plt.semilogx(w_ntf, h_ntf_db)
# plt.semilogx(f, Pxx)
# plt.ylim([-300, 100])
# plt.xlabel('frequency [Hz]')
# plt.ylabel('PSD [dB/Hz]')
# # plt.legend(["NTF freq response","error PSD"])
# plt.grid()
# plt.show()

# %%
