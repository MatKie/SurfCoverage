import numpy as np
from mkutils import save_to_file, create_fig
from Frumkin import FrumkinIsotherm

# +
# Curve fit
from scipy.optimize import curve_fit, leastsq
from scipy.optimize import root_scalar

experimental = np.loadtxt('Experimental.csv', delimiter=',')
pub_tens = np.loadtxt('Frumkin_tens.csv', delimiter=',')

def root_wrapper(x, Isotherm, c0):
    return Isotherm.concentration(x) - c0

def curve_fit_wrapper(concentrations, omega_0, b, alpha_f, epsilon):
    frumkin_params = {
        "omega_0": omega_0,
        "b": b,
        "alpha_f": alpha_f,
        "epsilon": epsilon,
        "zero_surf": 51.3,
    }
    tensions = []
    Isotherm = FrumkinIsotherm(**frumkin_params)
    gamma = 0.5
    concentration = concentrations[0]
    
    solution = root_scalar(root_wrapper, args=(Isotherm, concentration), x0 = gamma, x1= gamma+0.1)
    gamma = solution.root
    tension = Isotherm.tension(gamma)
    tensions.append(*tension)
    gammas = []
    gammas.append(gamma)
    for concentration0, concentration1 in zip(concentrations, concentrations[1:]):
        solution = root_scalar(root_wrapper, args=(Isotherm, concentration1), x0 = gamma, x1= gamma-0.01)
        if not solution.converged:
            for concentration in np.linspace(concentration0, concentration1, 1000):
                solution = root_scalar(root_wrapper, args=(Isotherm, concentration), x0 = gamma, x1= gamma-0.01)
                if not solution.converged:
                    tensions.append(1e4)
                    break
        else:
            gamma = solution.root
            tension = Isotherm.tension(gamma)
            gammas.append(gamma)
            tensions.append(*tension)
        
    return np.asarray(tensions)

def residuals(p, x, y):
    omega_min = (1.0 - p[3] * 51.3)
    if omega_min <= 1e-5:
        penalization = np.ones_like(x)*1e4
    else:
        penalization = np.zeros_like(x)
    return y - curve_fit_wrapper(x, p[0],p[1],p[2],p[3]) - penalization


#p0=[5.5e5, 0.532, -3.4, 0.009]
#p0, pcov = leastsq(func=residuals, x0=(pstart), args=(pub_tens[:-7, 0], pub_tens[:-7, 0]))

pstart = [5.5e+05,  532000, -3.4e+00,  0.009]
#pstart = [1e5, 1, -1, 0]
bounds = ([5.5e5-1e4, 0, -3.4-0.5, 0.009-0.01], [5.5e5+1e4, 1e6, -3.4+0.5, 0.009+0.005])
popt, pcov = curve_fit(curve_fit_wrapper, pub_tens[:-6, 0], pub_tens[:-6, 1], p0=pstart, bounds=bounds)
bounds = ([5.5e5-1e-8, 1e3, -3.4-1e-8, 0.009-1e-8], [5.5e5, 1e6, -3.4, 0.009])
popt_all, pcov = curve_fit(curve_fit_wrapper, pub_tens[:-6, 0], pub_tens[:-6, 1], p0=pstart, bounds=bounds)

#popt = [ 5.52403210e+05, 3.63788310e+01, -2.65209118e+00, 1.94931773e-02]

# +
pstart = [1.8e5, 0.53, -3.4e0, 0.0009]
bounds = ([1.8e5-1e-8, 0, -10, 0.0009], [1.8e5, 1e6, 10, 0.1])

p_bui, _ = curve_fit(curve_fit_wrapper, pub_tens[:-6, 0], pub_tens[:-6, 1], p0=pstart, bounds=bounds)
print(p_bui)

# +
print(popt,pcov)
print(popt, popt_all)

#popt = [ 5.39356572e+05,  3.00092656e+02, -2.98399513e+00,  9.95576496e-03]
curve_params = {
    "omega_0": popt[0],
    "b": popt[1],
    "alpha_f": popt[2],
    "epsilon": popt[3],
    "zero_surf": 51.3,
}
CurveFitIsotherm = FrumkinIsotherm(**curve_params)

#p0 = [ 5.62145355e+05,  3.97160931e+02, -2.72133622e+00,  9.74075322e-03]
#print(p0)
curve_all_params = {
    "omega_0": popt_all[0],
    "b": popt_all[1],
    "alpha_f": popt_all[2],
    "epsilon": popt_all[3],
    "zero_surf": 51.3,
}

CurveFitIsothermAll = FrumkinIsotherm(**curve_all_params)
#[ 5.62145355e+05  3.97160931e+02 -2.72133622e+00  9.74075322e-03]
# [ 5.39356572e+05  3.00092656e+02 -2.98399513e+00  9.95576496e-03] (to all)
llamas_params = {
    "omega_0": 5.5e5,
    "b": 5.32e-5,
    "alpha_f": -3.4,
    "epsilon": 0.009,
    "zero_surf": 51.3,
}
LlamasIsotherm = FrumkinIsotherm(**llamas_params)

# +
fig, ax = create_fig(1, 1)
ax = ax[0]
axt = ax.twinx()

fig1, ax1 = create_fig(1,1)
ax1 = ax1[0]
ax1t = ax1.twinx()

for Isotherm, ls, color in zip([CurveFitIsotherm, CurveFitIsothermAll, LlamasIsotherm], ['-', '--', 'dashdot'], ['C0', 'C1', 'C2']):
    gamma_max = Isotherm.get_gamma_zero_surftens()
    print(gamma_max)
    print(Isotherm.get_omega_zero_surftens())
    dx = 0.02
    gamma = np.linspace(0.00001 * gamma_max, (1 - dx) * gamma_max, num=5000)

    conc, tens = Isotherm.concentration(gamma), Isotherm.tension(gamma)
    ax1.plot(gamma, tens)
    ax1t.plot(gamma, conc)
    
    cap_omega = np.asarray([Isotherm.iterate_omega(gamma_i) for gamma_i in gamma])
    gamma = np.asarray([Isotherm.transform_gamma(gamma_i) for gamma_i in gamma])
    cap_omega = cap_omega * gamma
    ax.plot(
        conc, tens, lw=2, ls=ls
    )
    axt.plot(conc, cap_omega, lw=2, ls=ls, color=color)


ax.plot(experimental[:, 0], experimental[:, 1], ls='', marker='o')
ax.plot(pub_tens[:, 0], pub_tens[:, 1], ls='--',)
ax.set_xscale("log")
ax.set_xlabel("Concentration$\,/\,mol\,l^-1$")
ax.set_ylabel("Surface Tension$\,/\,mN\,m^-1$")
ax1.set_xlabel("$\Gamma\,/\,nm^2$")
ax1.set_ylabel("Surface Tension$\,/\,mN\,m^-1$")
axt.set_ylabel("$\Gamma \omega\,/\,-$")
ax1t.set_ylabel("Concentration$\,/\,mol\,l^-1$")
# +
gamma_max = Isotherm.get_gamma_zero_surftens()
dx = 0.1
gamma = np.linspace(0.00001 * gamma_max, (1 - dx) * gamma_max, num=1000)





fig, ax = create_fig(1, 1)
ax = ax[0]
omegas = [Isotherm.iterate_omega(gamma_i) for gamma_i in gamma]
ax.plot(Isotherm.concentration(gamma), np.asarray(omegas)*Isotherm.transform_gamma(gamma))
ax.set_xscale("log")
ax.set_ylabel("$\Gamma \omega\,/\,-$")
ax.set_xlabel("Concentration$\,/\,mol\,l^-1$")

fig, ax = create_fig(1, 1)
ax = ax[0]
omegas = [Isotherm.iterate_omega(gamma_i) for gamma_i in gamma]
#ax.plot(Isotherm.concentration(gamma), Isotherm.transform_gamma(gamma))
ax.plot(Isotherm.concentration(gamma), gamma)

ax.set_xscale("log")
ax.set_ylabel("$\Gamma,/\,-$")
ax.set_xlabel("Concentration$\,/\,mol\,l^-1$")

fig, ax = create_fig(1, 1)
ax = ax[0]
ax.plot(
    Isotherm.concentration(gamma), Isotherm.tension(gamma)
)
ax.plot(experimental[:, 0], experimental[:, 1], ls='', marker='o')
ax.plot(pub_tens[:, 0], pub_tens[:, 1], ls='--',)
ax.set_xscale("log")
ax.set_xlabel("Concentration$\,/\,mol\,l^-1$")
ax.set_ylabel("Surface Tension$\,/\,mN\,m^-1$")
    
if False:


    fig, ax = create_fig(1, 1)
    ax = ax[0]
    ax.plot(gamma, frumkin_conc(gamma, **frumkin_params))
    ax.set_yscale("log")
    ax.set_xlabel("$\Gamma\,/\,nm^2$")
    ax.set_ylabel("Concentration$\,/\,mol\,l^-1$")
# -





