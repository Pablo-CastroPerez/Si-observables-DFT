import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy.optimize import curve_fit

# Read data from the output file
lattice_parameters = []
total_energies = []

with open('si_energies', 'r') as file:
    # Skip the header
    next(file)
    for line in file:
        # Split each line and extract the lattice parameter and total energy
        parts = line.strip().split(',')
        if "No result found" not in line:
            parameter = float(parts[0].split('=')[1].strip().split(' ')[0])*10.26
            energy = float(parts[1].split(':')[1].strip().split(' ')[0])
            lattice_parameters.append(parameter)
            total_energies.append(energy)

lattice_parameters = np.array(lattice_parameters)
lattice_parameters = lattice_parameters*0.52918 
total_energies = np.array(total_energies) 

volumes_atom = lattice_parameters**3
volumes_atom = volumes_atom/8
total_energies_atom = total_energies/2
total_energies_atom = total_energies_atom * 13.605698 #eV

#EQUATIONS OF STATE

def murnaghan_eos(V, V0, B0, B0_prime, E0):
    eta = (V / V0)**(1/3)
    term1 = (B0 * V0 / B0_prime) * eta**3 * ((eta**(-3 * B0_prime)) / (B0_prime - 1) + 1)
    term2 = -(B0 * V0) / (B0_prime - 1)
    return E0 + term1 + term2

initial_guess_mur = [20, 100 , 4, -130]
popt_mur, pcov_mur = curve_fit(murnaghan_eos, volumes_atom, total_energies_atom, p0=initial_guess_mur, maxfev=20000)

# Extract the fitted parameters
V0_mur, B0_mur, B0_prime_mur, E0_mur = popt_mur

# Generate the fitted curve
volumes_fit_mur = np.linspace(min(volumes_atom), max(volumes_atom), 201)
energies_fit_mur = murnaghan_eos(volumes_fit_mur, V0_mur, B0_mur, B0_prime_mur, E0_mur)

residuals_mur = total_energies_atom - energies_fit_mur
ss_res_mur = np.sum((total_energies_atom - energies_fit_mur)**2)

# Calculate the total sum of squares (SS_tot)
ss_tot_mur = np.sum((total_energies_atom - np.mean(total_energies_atom))**2)

# Calculate R^2
r_squared_mur = 1 - (ss_res_mur / ss_tot_mur)

print(f"R^2 = {r_squared_mur}")

def birch_murnaghan_eos(V, V0, B0, B0_prime, E0):
    eta = (V / V0)**(1/3)
    term1 = (eta**-2 - 1)**3 * B0_prime
    term2 = (eta**-2 - 1)**2 * (6 - 4 * eta**-2)
    return E0 + (9 / 16) * B0 * V0 * (term1 + term2)


initial_guess_bm = [20, 100 , 4, -130]
popt_bm, pcov_bm = curve_fit(birch_murnaghan_eos, volumes_atom, total_energies_atom, p0=initial_guess_mur, maxfev=20000)

# Extract the fitted parameters
V0_bm, B0_bm, B0_prime_bm, E0_bm = popt_bm

# Generate the fitted curve
volumes_fit_bm = np.linspace(min(volumes_atom), max(volumes_atom), 201)
energies_fit_bm = murnaghan_eos(volumes_fit_bm, V0_bm, B0_bm, B0_prime_bm, E0_bm)

residuals_bm = (total_energies_atom - energies_fit_bm)
ss_res_bm = np.sum((total_energies_atom - energies_fit_bm)**2)

# Calculate the total sum of squares (SS_tot)
ss_tot_bm = np.sum((total_energies_atom - np.mean(total_energies_atom))**2)

# Calculate R^2
r_squared_bm = 1 - (ss_res_bm / ss_tot_bm)

print(f"R^2 = {r_squared_bm}")

def vinet_eos(V, V0, B0, B0_prime, E0):
    eta = (V / V0)**(1/3)
    term1 = 4 * B0 * V0 / ((B0_prime - 1) ** 2)
    term2 = 2 * B0 * V0 / ((B0_prime - 1) ** 2)
    exponent = np.exp((3/2) * (B0_prime - 1) * (1 - eta))
    return E0 + term1 + term2 * exponent * (3 * (B0_prime - 1) * (1 - eta) - 2)

initial_guess_vinet = [20, 100 , 4, -130]
popt_vinet, pcov_vinet = curve_fit(vinet_eos, volumes_atom, total_energies_atom, p0=initial_guess_mur, maxfev=20000)

# Extract the fitted parameters
V0_vinet, B0_vinet, B0_prime_vinet, E0_vinet = popt_vinet

# Generate the fitted curve
volumes_fit_vinet = np.linspace(min(volumes_atom), max(volumes_atom), 201)
energies_fit_vinet = vinet_eos(volumes_fit_vinet, V0_vinet, B0_vinet, B0_prime_vinet, E0_vinet)

residuals_vinet = total_energies_atom - energies_fit_vinet
ss_res_vinet = np.sum((total_energies_atom - energies_fit_vinet)**2)

# Calculate the total sum of squares (SS_tot)
ss_tot_vinet = np.sum((total_energies_atom - np.mean(total_energies_atom))**2)

# Calculate R^2
r_squared_vinet = 1 - (ss_res_vinet / ss_tot_vinet)
print(f"R^2 = {r_squared_vinet}")

errors_mur = np.sqrt(np.diag(pcov_mur))
errors_bm = np.sqrt(np.diag(pcov_bm))
errors_vinet = np.sqrt(np.diag(pcov_vinet))


# Print results for each EOS
print("Murnaghan EOS Parameters and Errors:")
for param, err in zip(popt_mur, errors_mur):
    print(f"{param:.5f} ± {err:.5f}")

print("\nBirch-Murnaghan EOS Parameters and Errors:")
for param, err in zip(popt_bm, errors_bm):
    print(f"{param:.5f} ± {err:.5f}")

print("\nVinet EOS Parameters and Errors:")
for param, err in zip(popt_vinet, errors_vinet):
    print(f"{param:.6f} ± {err:.6f}")



#OBSERVABLES

V0_mur = 20.45134 * 1e-30
B0_Pa_mur = 88.08* 1e9
mass_Si_kg = 28.0855 * 1.66054e-27 
rho_mur = ( mass_Si_kg) / V0_mur
c_l_mur = np.sqrt(B0_Pa_mur / rho_mur)
print(c_l_mur)

V0_bm =  20.44764* 1e-30
B0_Pa_bm = 88.47* 1e9
mass_Si_kg = 28.0855 * 1.66054e-27 
rho_bm = ( mass_Si_kg) / V0_bm
c_l_bm = np.sqrt(B0_Pa_bm / rho_bm)
print(c_l_bm)

V0_vinet =  20.445903* 1e-30
B0_Pa_vinet = 88.751* 1e9
mass_Si_kg = 28.0855 * 1.66054e-27 
rho_vinet = ( mass_Si_kg) / V0_vinet
c_l_vinet = np.sqrt(B0_Pa_vinet / rho_vinet)
print(c_l_vinet)

sigma_v_L_mur = (1 / 2) * c_l_bm * np.sqrt((0.004 /88.751 ) ** 2 + ( 0.00004 / 20.445903) ** 2)

import numpy as np
from scipy.constants import hbar, k as kB, pi

def debye_temperature(n, cL):
    """Calculate the Debye temperature for silicon."""
    # Compute the Debye wave vector kD (corrected formula)
    kD = (6 * pi**2 * n)**(1/3)

    # Compute the Debye frequency
    omegaD = cL * kD

    # Compute the Debye temperature
    thetaD = (hbar * omegaD) / kB

    return thetaD

print(debye_temperature(0.04890955415374*1e30, 6237.693409983))
print(kB)

def debye_temperature_error(thetaD, n, cL, delta_n, delta_cL):
 
   # Relative errors
    rel_error_n = delta_n / n
    rel_error_cL = delta_cL / cL

    # Total relative error in Theta_D
    rel_error_thetaD = rel_error_cL + (1/3) * rel_error_n

    # Absolute error in Theta_D
    delta_thetaD = thetaD * rel_error_thetaD

    return delta_thetaD
print(debye_temperature_error(679.1491109158055,0.0489095541*1e30, 6237.693409983635, 0.0000001 *1e30, 0.1))


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker  # Import ticker for fine control

# Create the main figure
plt.figure(figsize=(8, 8))
skip = 3
# Main plot: Data and fit
ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=3)  # Share x-axis with residuals
ax1.scatter(volumes_atom[::skip], total_energies_atom[::skip], color='black', label='Data', marker='+', s = 35)
ax1.plot(volumes_fit_mur[::skip], energies_fit_mur[::skip], color='r', linewidth = 2, linestyle = 'dotted',label='Murnaghan EOS', alpha = 0.6)
ax1.plot(volumes_fit_bm[::skip], energies_fit_bm[::skip], color='orange', linewidth = 2, linestyle = 'dashed', label='Birch-Murnaghan EOS', alpha = 0.6)
ax1.plot(volumes_fit_vinet[::skip], energies_fit_vinet[::skip], color='b',linewidth = 1, linestyle = 'solid',label='Vinet EOS',alpha = 0.6)
ax1.set_ylabel('Total Energy (eV)')  # Keep y-label
ax1.legend()
ax1.tick_params(axis='x', labelbottom=False)  # Remove x-axis labels
ax1.set_ylim(-130.580, -130.355)
ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

# Residual plot (inset or subplot)
ax2 = plt.subplot2grid((4, 1), (3, 0), rowspan=1, sharex=ax1)
ax2.scatter(volumes_atom[::skip], residuals_mur[::skip]*1000, color='r', label='Murhaghan Residuals', marker='d', s=10, alpha = 0.5)
ax2.scatter(volumes_atom[::skip], residuals_bm[::skip]*1000, color='orange', label='Birch-Murnaghan Residuals', marker='*', s=10,alpha = 0.6)
ax2.scatter(volumes_atom[::skip], residuals_vinet[::skip]*1000, color='b', label='Vinet Residuals', marker='s', s=7,alpha = 0.35)
ax2.axhline(0, color='k', linestyle='--', linewidth=1)  # Add a horizontal line at y=0
ax2.set_xlabel('Volume (Å$^3$)')
ax2.set_ylabel('Residuals (meV)', labelpad = 40)
ax2.legend()
ax2.set_ylim(-4.5, 9)
ax2.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))

# Reduce space between subplots to bring residuals closer to x-axis
plt.subplots_adjust(hspace=0)
plt.tight_layout()
plt.show()