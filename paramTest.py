#%%
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft

CLIGHT = 299792458
EPS0 = 8.85418782e-12
CHARGE_E = 1.602176634e-19
MASS_E = 9.10938356e-31

#%%
num_z = 8000
dz = 5.0e-9
dom_z = num_z * dz

print('num_z = {}'.format(num_z))
print('dz = {}'.format(dz))
print('dom_z = {}'.format(dom_z))

z1 = 10e-6
z2 = 30e-6
print('index of z1 = {}'.format(z1/dz))
print('index of z2 = {}'.format(z2/dz))

#%%
num_out = 22 - 1 
d_out = 1500
num_t = num_out * d_out
dt = dz / (2 * CLIGHT)
dom_t = num_t * dt

print('num_t = {}'.format(num_t))
print('dt = {}'.format(dt))
print('dom_t = {}'.format(dom_t))
print('dom_ct = {}'.format(CLIGHT * dom_t))

#%%
n0 = 1.5
domT_needed = dom_z * n0 / CLIGHT
numT_needed = domT_needed / dt
print('needed dom_t = {}'.format(domT_needed))
print('needed num_t = {}'.format(numT_needed))

#%%
num_freq = 1000
max_omeg = 2 * np.pi / ((num_t - 1) *dt) * num_freq
delta_om = 2 * np.pi / ((num_t - 1) *dt)

print('max_omeg = {:.3e}'.format(max_omeg))
print('delta_omeg = {:.3e}'.format(delta_om))

#%% Use the central frequency for further analysis
lamb0 = 1.9e-6
omeg0 = 2 * np.pi * CLIGHT / lamb0 
maxHarm = int(max_omeg / omeg0)

print('lamb0 = {}'.format(lamb0))
print('omeg0 = {}'.format(omeg0))
print('max harmonic = {:d}'.format(maxHarm))
print('time points per fund. wavelength = {:d}'.format(int((lamb0/CLIGHT)/dt)))
print('grid points per fund. wavelength = {:d}'.format(int(lamb0/dz)))

# %%
numE = 1.0e25
omegPlasma = np.sqrt(CHARGE_E**2*numE/(EPS0*MASS_E))

print("For {:} electrons, the plasma frequency is: {:.3e} 1/s".format(numE, omegPlasma))
print("Wavelength of {:.2e} [m] gives frequecy {:.2e} 1/s".format(lamb0, omeg0))
print("Approximate plasma skin depth: {:.2e}".format(CLIGHT/omegPlasma))
print("Critial Power: {:.3e} GW".format(17*(omeg0/omegPlasma)**2))
print("Frequecy (not angular) of central wavelength: {:.2e} Hz".format(omeg0/(2*np.pi)))

