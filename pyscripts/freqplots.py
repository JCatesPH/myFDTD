#%%
#import os, sys, subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#%%
mpl.rcParams['font.family'] = 'Tahoma'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
stdfigsize = (6.66, 5)

NOUT = 24
DOUT = 750 

#%% Read in output
print('Reading in grid.')
df = pd.read_csv('data/X.csv', header=None)
Z = df.values[:,0]
# %%
print("Reading in frequencies.")
df = pd.read_csv('data/omeg.csv', header=None)
omeg = df.values[:,0]

#%%
print('Reading in frequency domain array.')
df = pd.read_csv("data/Ew_re_0.csv", header=None, delimiter=',')
Ew_re = df.values
df = pd.read_csv("data/Ew_im_0.csv", header=None, delimiter=',')
Ew_im = df.values
Ew = Ew_re + 1.0j * Ew_im

#%%
print('Plotting frequency domain array.')
fig1 = plt.figure(figsize=(7,6), dpi=200)
ax1 = fig1.add_subplot(111)
C = ax1.contourf(omeg, Z, np.abs(Ew)**2)

ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.set_xlabel(r'$\omega$')
ax1.set_ylabel(r'$z$')
ax1.set_title(r'Intensity, $|E_z(\omega)|^2$')
fig1.colorbar(C)
plt.savefig('figs/Ew_0.png')

# %%
fig2 = plt.figure(figsize=(6,6), dpi=200)
ax2 = fig2.add_subplot(111)
line2, = ax2.plot(omeg, np.abs(Ew[int(Ew.shape[0]/2),:])**2)

ax2.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax2.set_xlabel(r'$\omega$')
ax2.set_title(r'Intensity spectrum, $|E_z(\omega)|^2$, at $z={z0:.2g}$'.format(z0=Z[int(Ew.shape[0]/2)]))
plt.savefig('figs/Ew_midpoint.png')


# %%
iom = 1
fig3 = plt.figure(figsize=(6,6), dpi=200)
ax3 = fig3.add_subplot(111)
line3, = ax3.plot(Z, np.abs(Ew[:,iom])**2)

ax3.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax3.set_xlabel(r'$z$')
ax3.set_title(r'Intensity spectrum, $|E_z(\omega)|^2$, at $\omega={omeg1:.2g}$'.format(omeg1=omeg[iom]))
plt.savefig('figs/Ew_omeg1.png')

#%%
for n in range(0, NOUT):
    T = n*DOUT

    df = pd.read_csv("data/Ew_re_{:d}.csv".format(T), header=None, delimiter=',')
    Ew_re = df.values
    df = pd.read_csv("data/Ew_im_{:d}.csv".format(T), header=None, delimiter=',')
    Ew_im = df.values
    Ew = Ew_re + 1.0j * Ew_im

    for coll in C.collections: # Remove the old data to replace
        plt.gca().collections.remove(coll) 

    C = ax1.contourf(omeg, Z, np.abs(Ew)**2)
    fig1.colorbar(C)
    fig1.canvas.draw()
    fig1.canvas.flush_events()
    fig1.savefig('figs/Ew_{:d}.png'.format(T))


    line2.set_ydata(np.abs(Ew[int(Ew.shape[0]/2),:])**2)
    fig2.canvas.draw()
    fig2.canvas.flush_events()
    fig2.savefig('figs/Ew_midpoint_{:d}.png'.format(T))

    line3.set_ydata(np.abs(Ew[:,iom])**2)
    fig3.canvas.draw()
    fig3.canvas.flush_events()
    fig3.savefig('figs/Ew_omeg1_{:d}.png'.format(T))

