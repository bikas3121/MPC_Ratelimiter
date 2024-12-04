
# %%
import numpy as np
from matplotlib  import pyplot as plt

beta_ideal = np.array([0, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01])
ENOB_ideal = np.array([8.39, 8.08, 7.72, 7.39, 6.96, 6.69, 6.01, 5.51, 5.37, 5.14])
f_sw_ideal = np.array([480, 396, 254, 218, 272, 310, 222.2, 236, 223, 201])


ENOB_inl = np.array([7.96, 7.73, 7.53, 7.28, 6.88, 6.68, 6.11, 5.55, 5.29, 5.13])
f_sw_inl = np.array([475, 379, 270, 278, 304, 303, 278, 249, 245, 229])

fig, axs = plt.subplots(2)
fig.suptitle('DAC-Ideal')
x = np.linspace(1,len(beta_ideal),len(beta_ideal))
axs[0].plot(x, ENOB_ideal, marker = "o")
axs[1].plot(x, f_sw_ideal, marker = "o")
axs[0].set_ylabel('ENOB')
axs[1].set_ylabel('$f_{sw}$ (kHz)' )
axs[1].set_xlabel(r'$\beta$')
axs[0].set_xticks(x, [])
axs[1].set_xticks(x)
axs[1].set_xticklabels(beta_ideal, rotation = 45)
axs[0].grid(axis = 'y')
axs[1].grid(axis = 'y')


fig, axs = plt.subplots(2)
fig.suptitle('DAC-INL')
x = np.linspace(1,len(beta_ideal),len(beta_ideal))
axs[0].plot(x, ENOB_inl, marker = "o")
axs[1].plot(x, f_sw_inl, marker = "o")
axs[0].set_ylabel('ENOB')
axs[1].set_ylabel('$f_{sw}$ (kHz)')
axs[1].set_xlabel(r'$\beta$')
axs[0].set_xticks(x, [])
axs[1].set_xticks(x)
axs[1].set_xticklabels(beta_ideal, rotation = 45)

axs[0].grid(axis = 'y')
axs[1].grid(axis = 'y')