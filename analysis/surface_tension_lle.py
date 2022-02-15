from mkutils import PlotGromacs, save_to_file, create_fig
import os, glob
import numpy as np
from Frumkin import FrumkinIsotherm

# +
Area = 8*8

def get_beads(coverage):
    return {'W2': 8976.-coverage, 'SO4V9': coverage, 'NA+': coverage, 'CM': 3.*coverage, 'CT': coverage }


# -

# Get energy files from all LLE simulations
files = glob.glob('../VLE_*')
files.sort(key=lambda x: int(x.split(os.sep)[1].split('_')[-1]))
files = files[:-3]
files = [os.path.join(file, 'prod_nvt', 'energies.out') for file in files]
files = [file for file in files if os.path.isfile(file)]
print([file.split(os.sep)[1].split('_')[-1] for file in files])


# + active=""
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
        print(coverage)
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
                # looks like avogadro * kb * 1000 so kJ/K/mol
                # Goes from kJ/mol to Kelvin
                correction /= 6.022*1380
        else:
            correction = 1.
        prop_averages.append(correction*data.get(prop).get('Average'))
        prop_errors.append(correction*data.get(prop).get('Error'))
        
    return coverage, prop_averages, prop_errors






# +
p0 = [ 5.50000000e+05,  9.53126622e+05, -3.40000000e+00,  9.00000000e-03]
llamas = [5.5e5, 0.532, -3.4, 0.009]

fig, ax = create_fig(1,1)
ax = ax[0]

for p in [p0, llamas]:
    frumkin_params = {
        "omega_0": p[0],
        "b": p[1],
        "alpha_f": p[2],
        "epsilon": p[3],
        "zero_surf": 51.3,
    }
    Isotherm = FrumkinIsotherm(**frumkin_params)
    gamma_max = Isotherm.get_gamma_zero_surftens()
    dx = 0.001
    gamma = np.linspace(0.001 * gamma_max, (1 - dx) * gamma_max, num=10000)


    coverage, surftens, surftens_err = get_data('#Surf*SurfTen', interaction_correction=False)

    ax.plot(gamma, Isotherm.tension(gamma), lw=2, ls='--', label='Frumkin Isotherm')  
    break
ax.errorbar(coverage, [prop_average/20. for prop_average in surftens], 
            [prop_error/20. for prop_error in surftens_err], marker='o', color='k', label='This Work', ls='')
#ax.plot(poly_cov(concentrations), poly_tens(concentrations), lw=2, ls='--', label='Two-State Model')

ax.set_xlabel('$\Gamma\,/\,nm^{-2}$')
ax.set_ylabel('Surface Tension$\,/\,mN\,m^{-1}$')
ax.legend()
save_to_file(os.path.join('LLE', 'Surface_tension'))


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
save_to_file(os.path.join('LLE', 'Coulomb'))

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
save_to_file(os.path.join('LLE', 'Water_interactions'))

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
save_to_file(os.path.join('LLE', 'SO4'))

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
save_to_file(os.path.join('LLE', 'CT'))
# -






