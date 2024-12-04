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
from MHOQ_RateLim import MHOQ_RLIM
from nsdcal import nsdcal, noise_shaping
from filter_parameters import filter_params

from process_sim_output import process_sim_output
from balreal import balreal
from welch_psd import welch_psd
from simulation_plots import bar_plot 
from figures_of_merit import TS_SINAD


# Generate test signal
def test_signal(SCALE, MAXAMP, FREQ, Rng,  OFFSET, t):
    Xcs = (SCALE/100)*MAXAMP*np.sin(2*np.pi*FREQ*t)   +Rng/2
    return Xcs 

# %% Choose MHOQ implemetatio method
class MHOQ_IMPLEMENTATION:
    BINARY = 1
    SCALED = 2
MHOQ_METHOD = MHOQ_IMPLEMENTATION.BINARY
# MHOQ_METHOD = MHOQ_IMPLEMENTATION.SCALED

# %% Quantiser configurations 
Qconfig = 2
Nb, Mq, Vmin, Vmax, Rng, Qstep, YQ, Qtype = quantiser_configurations(Qconfig, MHOQ_METHOD)

# %% Sampling frequency and rate
Fs = 1e6
Ts = 1/Fs

# %% Generate time vector
Xcs_SCALE = 100
Xcs_FREQ = 9999
match 2:
    case 1:  # specify duration as number of samples and find number of periods
        Nts = 1e5  # no. of time samples
        Np = np.ceil(Xcs_FREQ*Ts*Nts).astype(int) # no. of periods for carrier

    case 2:  # specify duration as number of periods of carrier
        Np = 50 # no. of periods for carrier
        
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
MLns_ME = YQns + INL_ME*Qstep   # Measured quantiser levels with measurement error

# %% Direct Quantization 
# C_DQ = quantise_signal(Xcs, Qstep, YQns, Qtype)
C_DQ = (np.floor(Xcs/Qstep +1/2)).astype(int)
match QMODEL:
    case 1:
        Xcs_DQ = generate_dac_output(C_DQ, YQns)
    case 2:
        Xcs_DQ = generate_dac_output(C_DQ, MLns_ME)

# %% NSD CAL
b_nsf = b_lpf - a_lpf
a_nsf = b_lpf
C_NSQ = noise_shaping(Nb, Xcs, b_nsf, a_nsf, Qstep, YQns, MLns, Vmin, QMODEL)  
match QMODEL:
    case 1:
        Xcs_NSQ = generate_dac_output(C_NSQ, YQns)
    case 2:
        Xcs_NSQ = generate_dac_output(C_NSQ, MLns_ME) 

# %% MHOQ Ratelimit
N_PRED = 1
Step = 0


# MHOQ_RL = MHOQ_RLIM(Nb, Qstep, QMODEL, A, B, C, D)
# C_MHOQ = MHOQ_RL.get_codes(N_PRED, Xcs, YQns, MLns, Step)
# match QMODEL:
#     case 1:
#         Xcs_MHOQ = generate_dac_output(C_MHOQ, YQns)
#     case 2:
#         Xcs_MHOQ = generate_dac_output(C_MHOQ, MLns_ME) 
# %%

# SWITCH_MIN = True
match QMODEL:
    case 1:
        QL = YQns.squeeze()
    case 2:
        QL = MLns.squeeze()
# Storage container for code
C_Store = []

# Loop length
len_MPC = Xcs.size - N_PRED

# len_MPC = 2
L =1e6
# State dimension
x_dim =  int(A.shape[0]) 

# Initial state
init_state = np.zeros(x_dim).reshape(-1,1)
steps = 2

# u_kminus1 = np.zeros((2**Nb, N_PRED))
u_kminus1_ind = (np.floor(Xcs[0]/Qstep +1/2)).astype(int) 
# u_kminus1_ind = 14 
# MPC loop
for j in tqdm.tqdm(range(len_MPC)):
# for j in range(1, len_MPC):
    m = gp.Model("MPC- INL")
    u = m.addMVar((2**Nb, N_PRED), vtype=GRB.BINARY, name= "u") # control variable
    x = m.addMVar((x_dim*(N_PRED+1),1), vtype= GRB.CONTINUOUS, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "x")  # State varible 
    # Add objective function
    Obj = 0

    # Set initial constraint
    m.addConstr(x[0:x_dim,:] == init_state)
    for i in range(N_PRED):
        k = x_dim * i
        st = x[k:k+x_dim]
        bin_con =  QL.reshape(1,-1) @ u[:,i].reshape(-1,1) 
        con = bin_con - Xcs[j+i]

        # Switching set 
        # Objective update
        e_t = C @ x[k:k+x_dim] + D * con
        Obj = Obj + e_t * e_t 

        # Constraints update
        f_value = A @ st + B * con
        st_next = x[k+x_dim:k+2*x_dim]
        m.addConstr(st_next == f_value)

        #% Limit search space
        ub = (u_kminus1_ind +  steps)
        if ub >= 2**Nb-1:
            ub = 2**Nb-1
        lb = (u_kminus1_ind -  steps)
        if lb <= 0:
            lb = 0        
        U1 = u[lb:ub+1,i].reshape(-1,1) 
        # Add constraint that sum of variables should be 1
        consU1 = gp.quicksum(U1[:,i]) == 1
        m.addConstr(consU1)

        # extract other variables that are not in the given bound
        mask = np.ones(u.size, dtype=bool)
        mask[np.arange(lb,ub+1,1)] = False
        U2 = u[mask]
        # Binary varialble constraint
        consU2 = U2== np.zeros([U2.shape[0],1])
        m.addConstr(consU2)

        # # Binary constraint
        # B1 = gp.quicksum(u[:,i]) == 1
        # m.addConstr(B1)
        # m.update
        
        # y_t = bin_con
        # y_t_minus1 = QL.reshape(1,-1) @ u_kminus1

        # lr1 = y_t - y_t_minus1
        # lr2 = -y_t + y_t_minus1
        # m.addConstr(lr1 <= L*Ts)
        # m.addConstr( lr2<= L* Ts )

    m.update
    # Set Gurobi objective
    m.setObjective(Obj, GRB.MINIMIZE)

    # 0 - Supress log output, 1- Print log outputs
    m.Params.LogToConsole = 0

    # Gurobi setting for precision  
    m.Params.IntFeasTol = 1e-9
    m.Params.IntegralityFocus = 1

    # Optimization 
    m.optimize()

    # Extract variable values 
    allvars = m.getVars()
    values = m.getAttr("X",allvars)
    values = np.array(values)

    # Variable dimension
    nr, nc = u.shape
    u_val = values[0:nr*nc]
    u_val = (np.reshape(u_val, (2**Nb, N_PRED))).astype(int)

    # Extract Code
    C_MPC = []
    for i in range(N_PRED):
        c1 = np.nonzero(u_val[:,i])[0][0]
        c1 = int(c1)
        C_MPC.append(c1)
    C_MPC = np.array(C_MPC)
    C_Store.append(C_MPC[0])

    U_opt = QL[C_MPC[0]] 
    # State prediction 
    con = U_opt - Xcs[j]
    x0_new =  A @ init_state + B * con
    # State update for subsequent prediction horizon 
    init_state = x0_new
    u_kminus1_ind = C_MPC[0]
    # j = j+1
C_MHOQ  = np.array(C_Store).reshape(1,-1)



match QMODEL:
    case 1:
        Xcs_MHOQ = generate_dac_output(C_MHOQ, YQns)
    case 2:
        Xcs_MHOQ = generate_dac_output(C_MHOQ, MLns_ME) 


# %% Trim vector lengths  
tm = t[:Xcs_MHOQ.shape[1]]
len_tm = len(tm)
TRANSOFF = np.floor(Npt*Fs/Xcs_FREQ).astype(int)  # remove transient effects from output

# %%  SINAD COmputation
SINAD_COMP_SEL = 1
plot_val = False

# FXcs_MHOQ, SINAD_MHOQ, ENOB_MHOQ= process_sim_output(tm, Xcs_MHOQ, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ")
# bar_plot(descr= "ENOB", MHOQ = ENOB_MHOQ)

FXcs_DQ, SINAD_DQ, ENOB_DQ= process_sim_output(tm, Xcs_DQ, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="NSQ")
FXcs_NSQ, SINAD_NSQ, ENOB_NSQ= process_sim_output(tm, Xcs_NSQ, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="NSQ")
FXcs_MHOQ1, SINAD_MHOQ1, ENOB_MHOQ1= process_sim_output(tm, Xcs_MHOQ, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ")
# FXcs_MHOQ2, SINAD_MHOQ2, ENOB_MHOQ2= process_sim_output(tm, Xcs_MHOQ2, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ")
# FXcs_MHOQ3, SINAD_MHOQ3, ENOB_MHOQ3= process_sim_output(tm, Xcs_MHOQ3, Fc, Fs, N_lp, TRANSOFF,SINAD_COMP_SEL,  plot=plot_val, descr="MHOQ")
# 
bar_plot(descr= "ENOB", DQ = ENOB_DQ,  NSQ = ENOB_NSQ, MHOQ1 = ENOB_MHOQ1)


 # %%
sl = 1000
fig, ax = plt.subplots()
ax.plot(t[10:sl], Xcs[10:sl])
ax.plot(t[10:sl], Xcs_MHOQ.squeeze()[10:sl])
# ax.plot(t[10:sl], Xcs_DQ.squeeze()[10:sl])
# ax.plot(t[10:sl], Xcs_NSQ.squeeze()[10:sl])
# ax.legend(['Ref','MHOQ',"DQ", "NSQ"])


# %%
