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
ETA0 = 376.730313668

#%%
lamb0 = 1.9e-6
omeg0 = 2 * np.pi * CLIGHT / lamb0 
n0 = 1.9

dz = 10.0e-9
z1 = 15e-6
z2 = 35e-6
zf = 50e-6

pml_len = lamb0 / 2 / dz
buff_len = lamb0 / dz
num_z = zf / dz + 2 * (pml_len + buff_len)
dom_z = num_z * dz
num_out = 50
d_out = 500
num_freq = 600

print('LZ >= {}'.format(pml_len))
print('BZ >= {}'.format(buff_len))
print('num_z = {}'.format(num_z))
print('dz = {}'.format(dz))
print('dom_z = {}'.format(dom_z))
print('index of z1 = {} + LZ + BZ'.format(z1/dz))
print('index of z2 = {} + LZ + BZ'.format(z2/dz))

#%%
num_t = (num_out - 1) * d_out
dt = dz / (2 * CLIGHT)
dom_t = num_t * dt

print('num_t = {}'.format(num_t))
print('dt = {}'.format(dt))
print('dom_t = {}'.format(dom_t))
print('dom_ct = {}'.format(CLIGHT * dom_t))

#%%
domT_needed = dom_z * n0 / CLIGHT
numT_needed = domT_needed / dt
print('needed dom_t = {}'.format(domT_needed))
print('needed num_t = {}'.format(numT_needed))

#%%
max_omeg = 2 * np.pi / ((num_t - 1) *dt) * num_freq
delta_om = 2 * np.pi / ((num_t - 1) *dt)

print('max_omeg = {:.3e}'.format(max_omeg))
print('delta_omeg = {:.3e}'.format(delta_om))

#%% Use the central frequency for further analysis
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

#################################################################################


# %%
