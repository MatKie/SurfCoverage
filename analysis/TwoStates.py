# -*- coding: utf-8 -*-
from scipy.interpolate import interp1d
import numpy as np
from mkutils import create_fig



# +
coverage = np.loadtxt('twostate_cov.csv', delimiter=',')
surftens = np.loadtxt('twostate_tens.csv', delimiter=',')
gamma = np.loadtxt('gamma_tot.csv', delimiter=',')

f_cov = interp1d(coverage[:, 0], coverage[:, 1], kind='cubic')
f_surf = interp1d(surftens[:, 0], surftens[:,1], kind='cubic')
f_gamma = interp1d(gamma[:, 0], gamma[:,1], kind='cubic')
f_gamma_t = lambda x: 6.022e5*f_gamma(x)

c_min = max([surftens[0,0], gamma[0,0]])
c_max = min([surftens[-1,0], gamma[-1,0]])
concentrations = np.geomspace(c_min, c_max)
# -



# +
fig, ax = create_fig(1,1)
ax = ax[0]
x = np.geomspace(coverage[0, 0], coverage[-1, 0], 100)
ax.plot(coverage[:, 0], coverage[:, 1], marker='x', ls='')
ax.plot(x, f_cov(x), ls = '-', lw=2)
ax.set_xlabel('Concentration / mol/l')
ax.set_ylabel('$\Gamma$ / mol/$m^2$')

ax.set_xscale('log')

fig, ax = create_fig(1,1)
ax = ax[0]
x = np.geomspace(surftens[0, 0], surftens[-1, 0], 100)
ax.plot(surftens[:, 0], surftens[:, 1], marker='x', ls='')
ax.plot(x, f_surf(x), ls = '-', lw=2)
ax.set_xlabel('Concentration / mol/l')
ax.set_ylabel('Interfacial Tension / mN/m')
ax.set_xscale('log')

fig, ax = create_fig(1,1)
ax = ax[0]
x = np.geomspace(gamma[0, 0], gamma[-1, 0], 100)
ax.plot(gamma[:, 0], gamma[:, 1], marker='x', ls='')
ax.plot(x, f_gamma(x), ls = '-', lw=2)
ax.set_xlabel('Concentration / mol/l')
ax.set_xlabel('$\Gamma$ / mol/$m^2$')

ax.set_xscale('log')

fig, ax = create_fig(1,1)
ax = ax[0]
ax.plot(f_gamma_t(concentrations), f_surf(concentrations))
ax.set_xlabel('$\Gamma$ / #/$nm^2$')
ax.set_ylabel('Interfacial Tension / mN/m')


# -


