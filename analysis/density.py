from mkutils import PlotGromacs, save_to_file, create_fig
import os, glob
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit, newton

# +
file = os.path.join('..', 'VLE_300', 'prod_nvt', 'density.xvg')

density = np.loadtxt(file, comments=['@', '#'])

labels = {1: 'Water', 2: 'NA+', 3: 'SO4', 4: 'CM', 5: 'CT', 6: 'CM_CT' }
# invert this dict.
indices = {v: k for k, v in labels.items()}

A = 8*8
box_size = 40.
def get_surfCov(file, A):
    return float(file.split(os.sep)[1].split('_')[-1])/A


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
file=files[-1]
print(file)


# +
def get_max(file, typ='SO4'):
    density = np.loadtxt(file, comments=['@', '#'])
    index = indices.get(typ)
    # Peaks are found from left to right.
    maxima, _ = find_peaks(density[:, index], distance=3, height=0.15*max(density[:, index]))
    
    #fig, ax = create_fig(1,1)
    #ax = ax[0]
    #ax.plot(density[:, 0], density[:, index], lw=2)
    #for i, maximum in enumerate(maxima):
    #    ax.axvline(density[maximum, 0], ls='--', color='C{:d}'.format(i+1))
    
    return [density[i, index] for i in maxima]



for typ in ['SO4', 'CM_CT']:
    fig, ax = create_fig(1,1)
    ax = ax[0]
    for file in files:
        surface_coverage = get_surfCov(file, A)
        maxima = get_max(file, typ=typ)

        if len(maxima) > 3:
            ax.plot(surface_coverage, np.average(maxima[1:3]), color='k', ls='', marker='o')
            ax.plot(surface_coverage, np.average(maxima[::3]), color='k', ls='', marker='o')
            ax.plot([surface_coverage]*2, [np.average(maxima[::3]), np.average(maxima[1:3])], color='k', ls='--')
        elif len(maxima) > 2:
            ax.plot(surface_coverage, min(maxima), color='k', ls='', marker='o')
            ax.plot(surface_coverage, (np.sum(maxima)-min(maxima))/2., color='k', ls='', marker='o')
            ax.plot([surface_coverage]*2, [min(maxima), (np.sum(maxima)-min(maxima))/2.], color='k', ls='--')

        else:
            ax.plot(surface_coverage, np.average(maxima), color='k', ls='', marker='o')
           
        ax.legend([typ])
        ax.set_xlabel('$\Gamma\,/\,nm^-2$')
        ax.set_ylabel('Peak $\\rho_N\,/\,nm^-2$')    
    save_to_file('Peak_Heights_{:s}'.format(typ))
# +
# Gaussian fits to CM_CT density and minimum value

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def get_lower_point(file, box_size=40, plot=False):
    density = np.loadtxt(file, comments=['@', '#'])
    index = indices.get('CM_CT')
    
    if plot:
        fig, ax = create_fig(1, 1)
        ax = ax[0]
        ax.plot(density[:, 0], density[:, index], lw=2)
    
    maxima, _ = find_peaks(density[:, index], distance=3, height=0.15*max(density[:, index]))
    points = []
    for maximum in maxima:
        # First get bounds for curve fitting
        compare_to = 0.01*density[maximum, index]
        array_length = len(density[:,1])
        for i, v in enumerate(density[maximum:, index]):
            if v < compare_to or maximum + i == array_length - 1:
                break
        right_start = maximum + i + 11

        for i, v in enumerate(density[maximum::-1, index]):
            if v < compare_to or maximum - i == 0:
                break
        left_start = maximum - i - 10
        
        # Lets fit a gaussian to this.. 
        p0 = [density[maximum, index], density[maximum, 0], 1] # a, x0, std
        z, _ = curve_fit(gauss_function, density[left_start:right_start, 0], density[left_start:right_start, index], p0)

        if plot:
            ax.plot(density[left_start:right_start, 0], gauss_function(density[left_start:right_start, 0], *z),
               lw=2, ls='--')
        
        # Get 10% value.. 
        cutoff = 0.1*z[0]
        # starting point needs to be slightly left for first peak, slightly right for second one..
        if density[maximum, 0] < box_size/2.:
            starting_point = z[1] - z[2]
        else:
            starting_point = z[1] + z[2]
        x = newton(lambda x: gauss_function(x, *z) - cutoff, starting_point)
        if plot:
            ax.axvline(x, color='k', lw=2, ls='--')
        points.append(x)
    return points
        
        
        
        

get_lower_point(file, plot=True, box_size=box_size)


# +
def get_water_density(file, box_size=40, slab_halfwidth=3, plot=False):
    density = np.loadtxt(file, comments=['@', '#'])
    width = slab_halfwidth
    upper = np.where(density[:, 0] > box_size/2. - width)
    lower = np.where(density[upper[0], 0] < box_size/2. + width)
    liq_maxima = upper[0][lower[0]]
    
    width = 5
    upper = np.where(density[:, 0] > box_size - width)
    lower = np.where(density[:, 0] < 0 + width)
    vap_maxima = np.concatenate((lower[0], upper[0]), axis=0)
    
    liq_density = np.average(density[liq_maxima[0], 1])
    vap_density = np.average(density[vap_maxima[0], 1])
    
    if plot:
        fig, ax = create_fig(1,1)
        ax = ax[0]
        
        ax.plot(density[:, 0], density[:, indices.get('Water')])
        ax.plot(density[liq_maxima, 0], density[liq_maxima, indices.get('Water')], lw=2, marker='x', ls='')
        ax.plot(density[vap_maxima, 0], density[vap_maxima, indices.get('Water')], lw=2, marker='x', ls='')
        
    return liq_density, vap_density

def water_density(x, x0, d, liq, vap, side='lhs', alpha=2.19722):
    av_dens = (liq+vap)/2.
    side_switch = {'lhs': -1, 'rhs': 1}.get(side)
    return av_dens - (av_dens*np.tanh(side_switch*alpha*(x-x0)/d))

def fit_water_density(file, box_size=40, side='lhs', plot=False):
    density = np.loadtxt(file, comments=['@', '#'])
    liq_density, vap_density = get_water_density(file, box_size)
    
    p0 = {'lhs': (box_size/3., 1.), 'rhs': (2*box_size/2., 1.)}.get(side)
    
    length = density[:, 0].shape[0]
    this_slice = {'lhs': slice(0,int(length/2)+1), 'rhs': slice(int(length/2.), length)}.get(side)
    
    p, _ = curve_fit(lambda x, x0, d: water_density(x, x0, d, liq_density, vap_density, side=side),
                    density[this_slice, 0], density[this_slice, indices.get('Water')], p0, bounds=([0, 1e-4], [box_size, 5]))
    if plot:
        fig, ax = create_fig(1, 1)
        ax = ax[0]
        ax.plot(density[this_slice, 0], density[this_slice, indices.get('Water')], marker='x')
        x = np.linspace(min(density[this_slice, 0]), max(density[this_slice, 0]), 1001)
        ax.plot(x, water_density(x, *p, liq_density, vap_density, side))
        ax.axvline(p[0], color='k')
        ax.axvline(p[0]-0.5*p[1], color='k', ls='--')
        ax.axvline(p[0]+0.5*p[1], color='k', ls='--')
        
    x = newton(lambda x: water_density(x, *p, liq_density, vap_density, side)-(0.9*liq_density), p[0])
    
    if plot:
        ax.axvline(x, color='C3', ls='--')
        
    side_switch = {'lhs': -1, 'rhs': 1}.get(side)
    upper_limit = -0.5*side_switch*p[1] + p[0]
    
    if abs(x-upper_limit) > 1e-2:
        raise RuntimeError ('newton and parameter from fit do not agree')
        
    d = p[1]
    return upper_limit, d
    
fit_water_density(file, plot=True)
fit_water_density(file, plot=True, side='rhs')

# +
fig, ax = create_fig(1, 1)
ax = ax[0]
surface_coverages, water_d, alkane_d = [], [], []
for file in files:
    surface_coverage = get_surfCov(file, A)
    surface_coverages.append(surface_coverage)
    water_left_x, water_left_d = fit_water_density(file, plot=False, box_size=box_size)
    water_right_x, water_right_d = fit_water_density(file, plot=False, side='rhs', box_size=box_size)
    alkane_left_x, alkane_right_x = get_lower_point(file, plot=False, box_size=box_size)
    
    left_d = water_left_x - alkane_left_x
    right_d = alkane_right_x - water_right_x
    
    water_d.append(water_left_d/2 + water_right_d/2.)
    alkane_d.append(left_d/2. + right_d/2.)
    
ax.plot(surface_coverages, water_d, marker='o', color='k', label='Water')
ax.plot(surface_coverages, alkane_d, marker='d', color='k', label='Alkane - Water')
    
ax.legend()
ax.set_xlabel('$\Gamma\,/\,nm^-2$')
ax.set_ylabel('Interface Thickness$\,/\,nm$')
save_to_file('Interface_Thickness')
# -




