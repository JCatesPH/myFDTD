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

#DATAPATH = 'DATA_boyle'
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

""" #%%
def makemp4(framerate=2, outname=None, rm=False):
    os.chdir('figs')
    bash1 = './makemp4.sh -f {0}'.format(framerate)
    if outname is not None:
        bash1 += ' -o {0}'.format(outname)
    if rm is True:
        bash1 += ' -r'
    subprocess.run(bash1.split(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True)

    os.chdir('..')

print(sys.argv)
if sys.argv[2] == 'true':
    rmopt = True
else:
    rmopt = False

"""
CLIGHT = 299792458
EPS0 = 8.85418782e-12
CHARGE_E = 1.602176634e-19
MASS_E = 9.10938356e-31

#%% Load in the simulation parameters
simParam_df = pd.read_table(DATAPATH + '/SimParameters.dat', 
    header=None, names=['Variable', 'Value', 'Description']
    )

# Remove trailing whitespace on variable names
simParam_df['Variable'] = simParam_df['Variable'].str.strip()

#%%
NOUT = int(simParam_df.loc[simParam_df['Variable'] == 'NOUT'].values[0,1])
DOUT = int(simParam_df.loc[simParam_df['Variable'] == 'DOUT'].values[0,1])
DELTA_Z = simParam_df.loc[simParam_df['Variable'] == 'DZ'].values[0,1]
DELTA_T = simParam_df.loc[simParam_df['Variable'] == 'DT'].values[0,1]
z1 = DELTA_Z * simParam_df.loc[simParam_df['Variable'] == 'iZ1'].values[0,1]
z2 = DELTA_Z * simParam_df.loc[simParam_df['Variable'] == 'iZ2'].values[0,1]

timearr = DELTA_T * np.arange(0, NOUT*DOUT, DOUT)

#%% Read in output
print('Reading in grid.')
df = pd.read_csv(DATAPATH + '/X.csv', header=None)
Z = df.values[:,0]


#%%
print('Reading in field array.')
df = pd.read_csv(DATAPATH + "/Ex.csv", header=None, delimiter=',')
Ez = df.values


# %%
print('Plotting field.')
plt.clf()
fig = plt.figure(figsize=(6,6), dpi=200)
ax = fig.add_subplot(111)
line1, = ax.plot(Z, Ez[2], 'b-')
ax.set_ylim(min(0.9*Ez.min(),1.1*Ez.min()), 1.1*Ez.max())
ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$E_z$')
""" text1 = ax.text(0.75, 0.9, 
    r'max$(E_z)=${:.3e}'.format(Ez[2].max()),
    horizontalalignment='center',
    verticalalignment='center',
    transform=ax.transAxes) """

ax.axvline(z1, color='k', linestyle='--')
ax.axvline(z2, color='k', linestyle='--')

#%%
for n in range(0, NOUT):
    T = n*DOUT

    line1.set_ydata(Ez[n])
    ax.set_title(r"$E_x(z)$ at step {t}".format(t=T))
    #text1.set_text(r'max$(E_z)=${:.3e}'.format(Ez[n].max()))
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.savefig(FIGPATH + 'output{:}.png'.format(n))


#print('Making animation.')
#makemp4(framerate=sys.argv[1], outname='test', rm=rmopt)



#%%
fig1 = plt.figure(figsize=(7,6), dpi=200)
ax1 = fig1.add_subplot(111)

#C = ax1.contourf(Z, timearr, np.abs(Ez)**2)
C = ax1.contourf(
        Z*1e6, 
        CLIGHT*timearr*1e6, 
        np.log10(np.abs(Ez)**2 + 1e-8),
        levels=50
        )

#ax1.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax1.set_xlabel(r'$z$ [$\mu$m]')
ax1.set_ylabel(r'$ct$ [$\mu$m]')
ax1.set_title(r'$|E(z,t)|^2$')
fig1.colorbar(C)

ax1.axvline(z1*1e6, color='w', linestyle='--')
ax1.axvline(z2*1e6, color='w', linestyle='--')

plt.savefig(FIGPATH + 'Ezt_contour.png')

#################################################################

#%%
print('Reading in carrier array.')
df = pd.read_csv(DATAPATH + "/rho.csv", header=None, delimiter=',')
rho = df.values

rho = np.where(rho < 1e-9, 1e-9, rho)
print('max(rho) = {:.2e}'.format(np.max(rho)))

#%%
fig1 = plt.figure(figsize=(7,6), dpi=200)
ax1 = fig1.add_subplot(111)

#C = ax1.contourf(Z, timearr, np.abs(Ez)**2)
C = ax1.contourf(
        Z*1e6, 
        CLIGHT*timearr*1e6, 
        np.log10(rho),
        levels=50
        )

#ax1.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax1.set_xlabel(r'$z$ [$\mu$m]')
ax1.set_ylabel(r'$ct$ [$\mu$m]')
ax1.set_title(r'$\rho(z,t)$')

ax1.axvline(z1*1e6, color='w', linestyle='--')
ax1.axvline(z2*1e6, color='w', linestyle='--')

fig1.colorbar(C)
plt.savefig(FIGPATH + 'rho_contour.png')
