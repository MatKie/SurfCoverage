from mkutils import PlotGromacs, save_to_file, create_fig
import os, glob
import numpy as np

# +
file = os.path.join('..', 'VLE_50', 'prod_nvt', 'density.xvg')

density = np.loadtxt(file, comments=['@', '#'])

labels = {1: 'Water', 2: 'NA+', 3: 'SO4', 4: 'CM', 5: 'CT', 6: 'CM_CT' }

# +
fig, ax = create_fig(1,1)
ax = ax[0]
ax2 = ax.twinx()

for i in [1,2,3,6]:
    ax.plot(density[:, 0], density[:, i], lw=2, label=labels.get(i), color='C{:d}'.format(i-1))

    
ax.set_xlim(10, 30)
ax.set_ylim(-0.5, max(density[:, -1])*1.1)   
ax.set_xlim()
ax.set_xlabel('z / nm')
ax.set_ylabel('$\\rho_N\ /\ 1\,/\,nm^3$')
ax.legend()

fig, ax = create_fig(1,1)
ax = ax[0]
for i in [4,5,6]:
    ax.plot(density[:, 0], density[:, i], lw=2, label=labels.get(i), color='C{:d}'.format(i-1))

    
ax.set_xlim(10, 30)
ax.set_ylim(-0.5, max(density[:, -1])*1.1)   
ax.set_xlim()
ax.set_xlabel('z / nm')
ax.set_ylabel('$\\rho_N\ /\ 1\,/\,nm^3$')
ax.legend()

# +
files = glob.glob('../VLE_*')
files = [os.path.join(file, 'prod_nvt', 'density.xvg') for file in files]
files = [file for file in files if os.path.isfile(file)]

files.sort(key=lambda x: int(x.split(os.sep)[1].split('_')[-1]))
files

# +
fig, ax = create_fig(1, 1)
ax = ax[0]

for file in files:
    density = np.loadtxt(file, comments=['@', '#'])
    ax.plot(density[:,0], density[:, 3], label=file.split(os.sep)[1].split('_')[-1])
    
ax.set_xlim(10, 30)
ax.legend()
# -


