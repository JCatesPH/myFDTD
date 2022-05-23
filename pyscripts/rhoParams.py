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

I0 = 1e12
taup = 20e-15
eDomT = 10.0 * taup
eNumT = 2500
t0 = eDomT / 2

sigma_k = 1e-23
k = 3
num_atoms = 2.1e26

timeArr = np.linspace(0, eDomT, eNumT)
rho = np.zeros_like(timeArr)
#%%
def eField(t):
    A0 = np.sqrt(2.0 * I0 / (CLIGHT * EPS0))
    return A0 * np.exp(-2.0 * np.log(2.0) * ((t-t0)/taup)**2) * np.cos(omeg0 * (t-t0))

#%%
for j in range(eNumT-1):
    intensity = eField(timeArr[j])**2 / (2.0*ETA0)
    drho = sigma_k * intensity**k * (num_atoms - rho[j])

    rho[j+1] = rho[j] + (eDomT / eNumT) * drho

    if (rho[j+1] > num_atoms):
        rho[j+1] = num_atoms

#%%
print('final rho = {:.3e}'.format(rho[-1]))
print('omega_p = {:.3e}'.format(np.sqrt(CHARGE_E**2*rho[-1]/(EPS0*MASS_E))))
print('omega ratio : {:.5f}'.format(np.sqrt(CHARGE_E**2*rho[-1]/(EPS0*MASS_E))/omeg0))

#%%
plt.plot(timeArr, eField(timeArr))
plt.show()

# %%
rho_filtered = np.where(rho < 1e-9, 1e-9, rho)
plt.semilogy(timeArr, rho_filtered)
