#%%
#import os, sys, subprocess
import sys
import glob
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#%%
mpl.rcParams['font.family'] = 'Tahoma'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
stdfigsize = (6.66, 5)

#%%
if len(sys.argv) > 1:
    DATAPATH = str(sys.argv[1])
else:
    DATAPATH = './data/'
if len(sys.argv) > 2:
    FIGPATH = str(sys.argv[2])
else:
    FIGPATH = './figs/'

# Set values if run in interactive mode (VSCode)
if hasattr(sys, 'ps1'):
    print("Interactive mode detected..")
    DATAPATH = '../data/'
    FIGPATH = '../figs/'

#%% Load in the simulation parameters
simParam_df = pd.read_table(DATAPATH + '/SimParameters.dat', 
    header=None, names=['Variable', 'Value', 'Description']
    )

# Remove trailing whitespace on variable names
simParam_df['Variable'] = simParam_df['Variable'].str.strip()

#%%
LZ = int(simParam_df.loc[simParam_df['Variable'] == 'LZ'].values[0,1])
BZ = int(simParam_df.loc[simParam_df['Variable'] == 'BZ'].values[0,1])
NOUT = int(simParam_df.loc[simParam_df['Variable'] == 'NOUT'].values[0,1])
DOUT = int(simParam_df.loc[simParam_df['Variable'] == 'DOUT'].values[0,1])
DELTA_Z = simParam_df.loc[simParam_df['Variable'] == 'DZ'].values[0,1]
DELTA_T = simParam_df.loc[simParam_df['Variable'] == 'DT'].values[0,1]
z1 = DELTA_Z * simParam_df.loc[simParam_df['Variable'] == 'iZ1'].values[0,1]
z2 = DELTA_Z * simParam_df.loc[simParam_df['Variable'] == 'iZ2'].values[0,1]

omeg0 = simParam_df.loc[simParam_df['Variable'] == 'omeg0'].values[0,1]

#%%
CLIGHT = 299792458
EPS0 = 8.85418782e-12
CHARGE_E = 1.602176634e-19
MASS_E = 9.10938356e-31

IFACTOR = EPS0*CLIGHT/2 / (NOUT * DOUT)


#%% Read in output
print('Reading in grid.')
df = pd.read_csv(DATAPATH + 'X.csv', header=None)
Z = df.values[:,0]
# %%
print("Reading in frequencies.")
df = pd.read_csv(DATAPATH + 'omeg.csv', header=None)
omeg = df.values[:,0]

#%%
print('Reading in frequency domain array.')
df = pd.read_csv(DATAPATH + "Ew_re.csv", header=None, delimiter=',')
Ew_re = df.values
df = pd.read_csv(DATAPATH + "Ew_im.csv", header=None, delimiter=',')
Ew_im = df.values
Ew = Ew_re + 1.0j * Ew_im

#%%
print('Plotting frequency domain array.')
fig1 = plt.figure(figsize=(7,6), dpi=200)
ax1 = fig1.add_subplot(111)

C = ax1.contourf(omeg/omeg0, Z, 
    np.log10(IFACTOR * np.abs(Ew)**2),
    levels=np.linspace(-1, 15, 50),
    extend='min'
    )
#C = ax1.contourf(omeg, Z, np.log10(np.abs(Ew)**2 + 1e-9))

ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.set_xlabel(r'$\omega/\omega_0$')
ax1.set_ylabel(r'$z$')
ax1.set_title(r'Intensity')
fig1.colorbar(C)
plt.savefig(FIGPATH + 'Ew.png')

# %%
fig2 = plt.figure(figsize=(6,6), dpi=200)
ax2 = fig2.add_subplot(111)
line2, = ax2.semilogy(omeg, IFACTOR * np.abs(Ew[int(Ew.shape[0]/2),:])**2)

#ax2.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
ax2.set_xlabel(r'$\omega$')
ax2.set_title(r'Intensity spectrum at $z={z0:.2g}$'.format(z0=Z[int(Ew.shape[0]/2)]))
plt.savefig(FIGPATH + 'Ew_midpoint.png')

# %%
plt.cla()
ax2.semilogy(omeg/omeg0, IFACTOR * np.abs(Ew[int(Ew.shape[0]) - (LZ+BZ),:])**2)
#line2.set_ydata(IFACTOR * np.abs(Ew[int(Ew.shape[0]) - 75,:])**2)
ax2.set_xlabel(r'$\omega/\omega_0$')
ax2.set_title(r"Transmitted")
fig2.canvas.draw()
fig2.canvas.flush_events()
fig2.savefig(FIGPATH + 'Ew_transmitted.png')

# %%
iom = 1
fig3 = plt.figure(figsize=(6,6), dpi=200)
ax3 = fig3.add_subplot(111)
line3, = ax3.plot(Z, IFACTOR * np.abs(Ew[:,iom])**2)

ax3.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax3.set_xlabel(r'$z$')
ax3.set_title(r'Intensity spectrum at $\omega={omeg1:.2g}$'.format(omeg1=omeg[iom]))
plt.savefig(FIGPATH + 'Ew_omeg1.png')

#%%
df = pd.read_csv(DATAPATH + "eRefl_re.csv", header=None, delimiter=',')
Erefl_re = df.values
df = pd.read_csv(DATAPATH + "eRefl_im.csv", header=None, delimiter=',')
Erefl_im = df.values
Erefl = np.array(Erefl_re + 1.0j * Erefl_im)

#%%
fig2 = plt.figure(figsize=(6,6), dpi=200)
ax2 = fig2.add_subplot(111)
ax2.semilogy(omeg/omeg0, IFACTOR * np.abs(Ew[int(Ew.shape[0]) - (LZ+BZ),:])**2, label='Transmitted')
ax2.semilogy(omeg/omeg0, IFACTOR * np.abs(Erefl.T)**2, label='Reflected')

#ax2.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
ax2.set_xlabel(r'$\omega/\omega_0$')
ax2.legend()
#ax2.set_title(r'')
plt.savefig(FIGPATH + 'Ew_both.png')

# %%
fig2 = plt.figure(figsize=(6,6), dpi=200)
ax2 = fig2.add_subplot(111)
ax2.semilogy(omeg/omeg0, IFACTOR * np.abs(Erefl.T)**2, label='Reflected')

#ax2.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
ax2.set_xlabel(r'$\omega/\omega_0$')
ax2.legend()
ax2.set_title(r'Reflected')
plt.savefig(FIGPATH + 'Ew_reflected.png')
