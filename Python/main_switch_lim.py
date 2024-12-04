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
import gurobipy as gp
from gurobipy import GRB


from static_dac_model import quantise_signal, generate_code, generate_dac_output
from  quantiser_configurations import quantiser_configurations, get_measured_levels

from MHOQ_FreqMIN import MHOQ_FMIN
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
    Xcs = (SCALE/100)*MAXAMP*np.cos(2*np.pi*FREQ*t)   +Rng/2
    return Xcs 

# %% Choose MHOQ implemetatio method
class MHOQ_IMPLEMENTATION:
    BINARY = 1
    SCALED = 2
MHOQ_METHOD = MHOQ_IMPLEMENTATION.BINARY
# MHOQ_METHOD = MHOQ_IMPLEMENTATION.SCALED

# %% Quantiser configurations 
Qconfig = 1
Nb, Mq, Vmin, Vmax, Rng, Qstep, YQ, Qtype = quantiser_configurations(Qconfig, MHOQ_METHOD)

# %% Sampling frequency and rate
Fs = 1e6
Ts = 1/Fs

# %% Generate time vector
Xcs_SCALE = 100
Xcs_FREQ = 5
match 2:
    case 1:  # specify duration as number of samples and find number of periods
        Nts = 1e5  # no. of time samples
        Np = np.ceil(Xcs_FREQ*Ts*Nts).astype(int) # no. of periods for carrier

    case 2:  # specify duration as number of periods of carrier
        Np = 1 # no. of periods for carrier
        
Npt = 1  # no. of carrier periods to use to account for transients
Np = Np + Npt

t_end = Np/Xcs_FREQ  # time vector duration
t = np.arange(0, t_end, Ts)  # time vector

# %% Generate carrier/test signal
SIGNAL_MAXAMP = Rng/2
SIGNAL_OFFSET = -Qstep/2  # try to center given quantiser type
Xcs = test_signal(Xcs_SCALE, SIGNAL_MAXAMP, Xcs_FREQ, Rng,  SIGNAL_OFFSET, t)

fig, ax = plt.subplots()
ax.plot(t,Xcs)
# %% Reconstruction filter
N_lp = 2        # Filter order
FcS = 1e5       # Filter cutoff frequency
b_r, a_r = signal.butter(N_lp, FcS/(Fs/2),"low")

# %% Choose filter type for filter parameters to use in MHOQ impelmentation
class RCF: # 
    Perception = 1          # Psychoacoustically optimal noise-shaping filter form Goodwin et. al 2003
    Butterworth = 2         # Butterworth filter
    Butterworth_H2_Hinif = 3 # ||NTF||_{2} <= 1.5           # Optimal filter obtained minimising the error variance. 
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
QMODEL = 2

# %% Quantiser levels : Ideal and Measured
YQns = YQ.squeeze()
INL  =  get_measured_levels(Qconfig).squeeze()
MLns = YQns + INL*Qstep # Measured quantiser levels without measurement error

# Add measurement error to INL
sd_inl = 0.01*statistics.stdev(INL)
Measr_err = np.round(np.random.uniform(-sd_inl,sd_inl,len(INL)),4)
INL_ME = INL + Measr_err
MLns_ME = YQns + INL_ME*Qstep   # Measured quantiser levels with measurement error



# %% 
# %% Direct Quantization 
# C_DQ = quantise_signal(Xcs, Qstep, YQns, Qtype)
C_DQ = (np.floor(Xcs/Qstep +1/2)).astype(int)
match QMODEL:
    case 1:
        Xcs_DQ = generate_dac_output(C_DQ, YQns)
    case 2:
        Xcs_DQ = generate_dac_output(C_DQ, MLns)
# %% NSD CAL
b_nsf = b_lpf - a_lpf
a_nsf = b_lpf
C_NSQ = noise_shaping(Nb, Xcs, b_nsf, a_nsf, Qstep, YQns, MLns_ME, Vmin, QMODEL)  
match QMODEL:
    case 1:
        Xcs_NSQ = generate_dac_output(C_NSQ, YQns)
    case 2:
        Xcs_NSQ = generate_dac_output(C_NSQ, MLns) 
 # %% MHOQ
beta = 0.0
MHOQ_FreqMIN= MHOQ_FMIN(Nb, Qstep, QMODEL, A, B, C, D, beta)

N1 = 1
C_MHOQ1= MHOQ_FreqMIN.get_codes(N1, Xcs, YQns, MLns)
match QMODEL:
    case 1:
        Xcs_MHOQ1 = generate_dac_output(C_MHOQ1, YQns)
    case 2:
        Xcs_MHOQ1 = generate_dac_output(C_MHOQ1, MLns) 

# %% switch counter
S_counter  = 0
for i in range(1, Xcs_MHOQ1.size):
    if Xcs_MHOQ1[0,i-1] != Xcs_MHOQ1[0,i]:
        S_counter += 1
# print("Total number of switches:",S_counter)


# %% Average switching frequency
M = len(Xcs)
f_sw = (S_counter/(M*Ts*1e3))
print("Average switching frequency:",f_sw)
# %% Trim vector lengths  
tm = t[:Xcs_MHOQ1.shape[1]]
len_tm = len(tm)
TRANSOFF = np.floor(Npt*Fs/Xcs_FREQ).astype(int)  # remove transient effects from output

# %%  SINAD COmputation
SINAD_COMP_SEL = 1
plot_val = False
FXcs_DQ, SINAD_DQ, ENOB_DQ= process_sim_output(tm, Xcs_DQ, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="NSQ")
FXcs_NSQ, SINAD_NSQ, ENOB_NSQ= process_sim_output(tm, Xcs_NSQ, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="NSQ")
FXcs_MHOQ1, SINAD_MHOQ1, ENOB_MHOQ1= process_sim_output(tm, Xcs_MHOQ1, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ")
bar_plot(descr= "ENOB", DQ = ENOB_DQ,  NSQ = ENOB_NSQ, MHOQ1 = ENOB_MHOQ1)
# %%

# FXcs_MHOQ = signal.lfilter(b_r,a_r, Xcs_MHOQ1)
sl = 39000
fig, ax = plt.subplots()
ax.plot(t[10:sl], Xcs[10:sl])
ax.plot(t[10:sl], Xcs_MHOQ1.squeeze()[10:sl])
# ax.plot(t[10:sl], FXcs_MHOQ.squeeze()[10:sl])
# %%
