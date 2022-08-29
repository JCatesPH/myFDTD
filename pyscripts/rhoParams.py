#%%
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from numpy.fft import fft, ifft

CLIGHT = 299792458
EPS0 = 8.85418782e-12
CHARGE_E = 1.602176634e-19
MASS_E = 9.10938356e-31
ETA0 = 376.730313668


intensityFactor = EPS0*CLIGHT/2 * 1e-4 # Second factor changes it to W/cm^2

mpl.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 2
stdfigsize = (6.66, 5)

#%%
lamb0 = 1.9e-6
omeg0 = 2 * np.pi * CLIGHT / lamb0 

I0 = 1e12
taup = 20e-15
eDomT = 10.0 * taup
eNumT = 16000
t0 = eDomT / 2

sigma_k = 1e-23
k = 3
num_atoms = 2.0e26

timeArr = np.linspace(0, eDomT, eNumT)
rho = np.zeros_like(timeArr)
rho_env = np.zeros_like(timeArr)

#%%
def eEnvelope(t):
    A0 = np.sqrt(2.0 * I0 / (CLIGHT * EPS0))
    return A0 * np.exp(-2.0 * np.log(2.0) * ((t-t0)/taup)**2)

def eField(t):
    return eEnvelope(t) * np.cos(omeg0 * (t-t0))

#%%
for j in range(eNumT-1):
    intensity = eField(timeArr[j])**2 / (2.0*ETA0)
    drho = sigma_k * intensity**k * (num_atoms - rho[j])

    rho[j+1] = rho[j] + (eDomT / eNumT) * drho

    if (rho[j+1] > num_atoms):
        rho[j+1] = num_atoms

    intensity = eEnvelope(timeArr[j])**2 / (2.0*ETA0)
    drho = sigma_k * intensity**k * (num_atoms - rho_env[j])

    rho_env[j+1] = rho_env[j] + (eDomT / eNumT) * drho

    if (rho_env[j+1] > num_atoms):
        rho_env[j+1] = num_atoms

#%%
print('final rho = {:.3e}'.format(rho[-1]))
print('omega_p = {:.3e}'.format(np.sqrt(CHARGE_E**2*rho[-1]/(EPS0*MASS_E))))
print('omega ratio : {:.5f}'.format(np.sqrt(CHARGE_E**2*rho[-1]/(EPS0*MASS_E))/omeg0))

#%%
fig = plt.figure(figsize=(6,6), dpi=300)
ax = fig.add_subplot(111)

ax.plot(timeArr, eField(timeArr))
ax.plot(timeArr, eEnvelope(timeArr))

ax.margins(x=0)

plt.show()

# %%
rho_filtered = np.where(rho < 1e-9, 1e-9, rho)
rho_env_filtered = np.where(rho_env < 1e-9, 1e-9, rho_env)

fig = plt.figure(figsize=(6,6), dpi=300)
ax = fig.add_subplot(111)

ax.semilogy(timeArr, rho_filtered, label='Analytic signal')
ax.semilogy(timeArr, rho_env_filtered, label='Envelope')

ax.legend()
ax.margins(x=0)

plt.show()

#%%
distance = rho - rho_env

fig = plt.figure(figsize=(6,6), dpi=300)
ax = fig.add_subplot(111)

ax.plot(timeArr, distance)

ax.margins(x=0)

plt.show()


#%%
""" 
double atomicDensityProfile(int index) {
	double z0 = idie2 * DZ;
	double n_max = 2.0e26;

    double z = index * DZ;
	return (n_max / M_PI) / (pow(z - z0, 2) + pow(n_max, 2));
} 
"""

def atomicDensity(z):
    z0 = 10e-6
    n_max = 2.1e26
    gamma = 0.75e-9

    return n_max * (gamma**2  / ((z-z0)**2 + gamma**2))

#%%
zArr = np.linspace(0.0, 20e-6, 1000)

plt.subplots(figsize=stdfigsize, dpi=400)

plt.semilogy(zArr, atomicDensity(zArr)*1e-6)
plt.ticklabel_format(axis='x', style='sci', scilimits=(-6,-6))
plt.xlabel(r'$z$ [m]')
plt.ylabel(r'$N_a$ [cm$^{-3}$]')

plt.savefig('../figs/atomicProfile.png',dpi=400, transparent=False, bbox_inches='tight')
# %%
