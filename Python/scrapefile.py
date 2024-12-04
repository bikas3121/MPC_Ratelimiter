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
from dac_slew import rate_limiter



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
Qconfig = 1
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
        Np = 5 # no. of periods for carrier
        
Npt = 1  # no. of carrier periods to use to account for transients
Np = Np + Npt

t_end = Np/Xcs_FREQ  # time vector duration
t = np.arange(0, t_end, Ts)  # time vector

# %% Generate carrier/test signal
SIGNAL_MAXAMP = Rng/2
SIGNAL_OFFSET = -Qstep/2  # try to center given quantiser type
Xcs = test_signal(Xcs_SCALE, SIGNAL_MAXAMP, Xcs_FREQ, Rng,  SIGNAL_OFFSET, t)

QMODEL = 1
YQns = YQ
# %% Direct Quantization 
C_DQ = (np.floor(Xcs/Qstep +1/2)).astype(int)
match QMODEL:
    case 1:
        Xcs_DQ = generate_dac_output(C_DQ, YQns)
    case 2:
        Xcs_DQ = generate_dac_output(C_DQ, MLns_ME)

# %% Rate limiting 
Xcs_rls = Xcs_DQ.squeeze()
Xcs_RL = np.zeros_like(Xcs_rls)

u0 = Xcs_rls[i]
for i in range(len(Xcs)-1):
    u1 = Xcs_rls[i+1]
    t0 = t[i]
    t1 = t[i+1]
    R = 1e3
    F = 1e4
    u0 = rate_limiter(u0, u1, t0,t1,R, F) 
    Xcs_RL[i] = u0

# %%
sl = 20
fig, ax = plt.subplots()
ax.plot(t[:sl],Xcs_RL[:sl])
ax.plot(t[:sl],Xcs_DQ.squeeze()[:sl])

# %%
