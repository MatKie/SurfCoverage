from mkutils import PlotGromacs, save_to_file, create_fig
import os, glob

Area = 8.*8. # nm

files = glob.glob('../VLE_*')
files = [os.path.join(file, 'prod_nvt', 'energies.out') for file in files]
files = [file for file in files if os.path.isfile(file)]
print(files)


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

def get_data(prop):
    coverage, prop_averages, prop_errors = [], [], []

    for file in files:
        data = PlotGromacs.get_gmx_stats(file)
        prop_averages.append(data.get(prop).get('Average'))
        prop_errors.append(data.get(prop).get('Error'))
        coverage.append(float(file.split(os.sep)[1].split('_')[-1])/Area)
    return coverage, prop_averages, prop_errors



# +
fig, ax = create_fig(1,1)
ax = ax[0]
coverage, surftens, surftens_err = get_data('#Surf*SurfTen')
ax.errorbar(coverage, [prop_average/20. for prop_average in surftens], 
            [prop_error/20. for prop_error in surftens_err], ls='', marker='o')

ax.set_xlabel('Headgroup per Area / 1/nm$^2$')
ax.set_ylabel('Surface Tension / Nm')

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
coverage, prop_averages, prop_errors = get_data(prop)
ax.errorbar(coverage, [prop_average for prop_average in prop_averages], 
            [prop_error for prop_error in prop_errors], ls='', marker='o', color='C3', label=prop)

ax.set_xlabel('Headgroup per Area / 1/nm$^2$')
ax.set_ylabel('Interactions strength / kJ/mol')
ax.legend()

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
# -






