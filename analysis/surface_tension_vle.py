from mkutils import PlotGromacs, save_to_file, create_fig
import os, glob
import numpy as np

# +
Area = 8*8

def get_beads(coverage):
    return {'W2': 8976.-coverage, 'SO4V9': coverage, 'NA+': coverage, 'CM': 3.*coverage, 'CT': coverage }


# -

files = glob.glob('../VLE_*')
files.sort(key=lambda x: int(x.split(os.sep)[1].split('_')[-1]))
#files = files[:-3]
files = [os.path.join(file, 'prod_vle', 'energies.out') for file in files]
files = [file for file in files if os.path.isfile(file)]
print([file.split(os.sep)[1].split('_')[-1] for file in files])


# + active=""
#
# dict_keys(['Bond', 'Angle', 'LJ (SR)', 'Coulomb (SR)', 'Coul. recip.', 'Potential', 'Kinetic En.', 'Total Energy', 
#            'Conserved En.', 'Temperature', 'Pressure', 'Pres-XX', 'Pres-YY', 'Pres-ZZ', '#Surf*SurfTen', 'LJ-SR:W2-W2', 
#            'LJ-SR:W2-CM', 'LJ-SR:W2-CT', 'LJ-SR:W2-NA+', 'Coul-SR:SO4V9-SO4V9', 'LJ-SR:SO4V9-SO4V9', 'LJ-SR:SO4V9-CM', 
#            'LJ-SR:SO4V9-CT', 'Coul-SR:SO4V9-NA+', 'LJ-SR:SO4V9-NA+', 'LJ-SR:CM-CM', 'LJ-SR:CM-CT', 'LJ-SR:CM-NA+', 
#            'LJ-SR:CT-CT', 'LJ-SR:CT-NA+', 'Coul-SR:NA+-NA+', 'LJ-SR:NA+-NA+'])
#            
# prop = '#Surf*SurfTen'
# prop = 'LJ-SR:SO4V9-NA+'
# prop = 'Pres-ZZ'
# -

def get_data(prop, interaction_correction=True, in_kelvin=False):
    coverage, prop_averages, prop_errors = [], [], []

    for file in files:
        data = PlotGromacs.get_gmx_stats(file)
        coverage.append(float(file.split(os.sep)[1].split('_')[-1])/Area)
        if interaction_correction:
            if prop == 'Coul. recip.':
                type_1, type_2 = 'NA+', 'SO4V9'
            else:
                interaction_string = prop.split(':')[-1]
                type_1, type_2 = interaction_string.split('-')
            beads = get_beads(coverage[-1])
            beads_1 = beads.get(type_1)
            beads_2 = beads.get(type_2)
            correction = 1./(beads_1*beads_2)
            if in_kelvin:
                correction /= 6.022*1380
        else:
            correction = 1.
        prop_averages.append(correction*data.get(prop).get('Average'))
        prop_errors.append(correction*data.get(prop).get('Error'))
        
    return coverage, prop_averages, prop_errors


# +
from scipy.interpolate import interp1d
gamma = np.loadtxt('posser_fig1_gamma.csv', delimiter=',')
surftens = np.loadtxt('posser_fig1_tens.csv', delimiter=',')

f_gamma = interp1d(gamma[:, 0], gamma[:, 1], kind='cubic')
f_surf = interp1d(surftens[:, 0], surftens[:,1], kind='cubic')
f_gamma_t = lambda x: 6.022e5*f_gamma(x)

c_min = max([surftens[0,0], gamma[0,0]])
c_max = min([surftens[-1,0], gamma[-1,0]])
concentrations = np.geomspace(c_min, c_max)
concentrations = surftens[1:-1, 0]

# +
gamma = np.loadtxt('posser_fig1_gamma.csv', delimiter=',')
surftens = np.loadtxt('posser_fig1_tens.csv', delimiter=',')

f_gamma = np.polyfit(gamma[:, 0], gamma[:, 1], deg=2)
f_gamma = np.poly1d(f_gamma)
f_surf = np.polyfit(surftens[:, 0], surftens[:,1], deg=2)
f_surf = np.poly1d(f_surf)
f_gamma_t = lambda x: 6.022e5*f_gamma(x)

c_min = min([surftens[0,0], gamma[0,0]])
c_max = max([surftens[-1,0], gamma[-1,0]])
concentrations = np.geomspace(c_min, c_max)
concentrations = surftens[:, 0]

# +
fig, ax = create_fig(2, 1, sharex=True)
ax1 = ax[-1]
ax = ax[0]


ax.plot(surftens[:, 0], surftens[:, 1], marker='o', ls='', color='k', fillstyle='none', markersize=8)
ax.plot(concentrations, f_surf(concentrations), lw=2)
ax1.plot(gamma[:, 0], gamma[:, 1], marker='o', ls='', color='k', fillstyle='none', markersize=8)
ax1.plot(concentrations, f_gamma(concentrations), lw=2)


# +
fig, ax = create_fig(1,1)
ax = ax[0]
coverage, surftens, surftens_err = get_data('#Surf*SurfTen', interaction_correction=False)
ax.errorbar(coverage, [prop_average/20. for prop_average in surftens], 
            [prop_error/20. for prop_error in surftens_err], ls='', marker='o', color='k')

ax.plot(f_gamma_t(concentrations), f_surf(concentrations), marker='', ls='--', lw=2)
ax.set_xlabel('$\Gamma\,/\,nm^{-2}$')
ax.set_ylabel('Surface Tension$\,/\,mN\,m^{-1}$')
save_to_file(os.path.join('VLE', 'Surface_tension'))

# +
fig, ax = create_fig(1,1, )
ax = ax[0]
prop = 'LJ-SR:SO4V9-NA+'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)


prop = "Coul-SR:SO4V9-NA+"
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', color='C1', label=prop)

prop = "Coul. recip."
coverage, prop_averages, prop_errors = get_data(prop, interaction_correction=True)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', color='C3', label=prop)

ax.set_xlabel('Headgroup per Area / 1/nm$^2$')
ax.set_ylabel('Interactions strength / kJ/mol')
ax.legend()
save_to_file(os.path.join('VLE', 'Coulomb'))

# +
fig, ax = create_fig(1,1, )
ax = ax[0]
prop = 'LJ-SR:W2-CM'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:W2-CT'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:W2-NA+'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:W2-SO4V9'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)




ax.set_xlabel('Headgroup per Area / 1/nm$^2$')
ax.set_ylabel('Interactions strength / kJ/mol')
ax.legend()
save_to_file(os.path.join('VLE', 'Water_interactions'))

# +
fig, ax = create_fig(1,1, )
ax = ax[0]
prop = 'LJ-SR:SO4V9-CM'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:SO4V9-CT'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:SO4V9-NA+'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:SO4V9-SO4V9'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)
prop = 'Coul-SR:SO4V9-SO4V9'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'Coul-SR:SO4V9-NA+'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)


ax.set_xlabel('Headgroup per Area / 1/nm$^2$')
ax.set_ylabel('Interactions strength / kJ/mol')
ax.legend()
save_to_file(os.path.join('VLE', 'SO4'))

# +
fig, ax = create_fig(1,1, )
ax = ax[0]
prop = 'LJ-SR:CM-CM'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:CM-CT'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:CM-NA+'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:SO4V9-CM'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:W2-CM'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)


ax.set_xlabel('Headgroup per Area / 1/nm$^2$')
ax.set_ylabel('Interactions strength / kJ/mol')
ax.legend()
save_to_file(os.path.join('LLE', 'CM'))

# +
fig, ax = create_fig(1,1, )
ax = ax[0]


prop = 'LJ-SR:CT-CT'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:CM-CT'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:CT-NA+'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:SO4V9-CT'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)

prop = 'LJ-SR:W2-CT'
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', label=prop)


ax.set_xlabel('Headgroup per Area / 1/nm$^2$')
ax.set_ylabel('Interactions strength / kJ/mol')
ax.legend()
save_to_file(os.path.join('VLE', 'CT'))
# -






