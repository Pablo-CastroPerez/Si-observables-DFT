import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score  
import matplotlib.ticker as ticker

# Initialize lists
ca_ratios = []
stress_xx = []
stress_yy = []
stress_zz = []

# Open and read the file
file_path = "c_over_a_stress_results.out"

with open(file_path, "r") as file:
    for line in file:
        line = line.strip()
        if "No result found" in line or "c/a Ratio" in line or not line:
            continue  # Skip missing or header lines

        values = [v.strip() for v in line.split(",")]

        if len(values) < 7 or any(v == "" for v in values[0:7]):
            continue  # Skip invalid lines

        ca_ratios.append(float(values[0]))
        stress_xx.append(float(values[4]))
        stress_yy.append(float(values[5]))
        stress_zz.append(float(values[6]))

# Convert to NumPy arrays
conversion_factor = 14710.5  # Ry/bohr³ to GPa

ca_ratios = np.array(ca_ratios)
strain = ca_ratios - 1
stress_xx = np.array(stress_xx)
stress_xx = - stress_xx
stress_xx = stress_xx*conversion_factor
stress_yy = np.array(stress_yy)

stress_yy = - stress_yy
stress_yy = stress_yy*conversion_factor

stress_zz = np.array(stress_zz)
stress_zz = - stress_zz
stress_zz = stress_zz*conversion_factor

# Define a linear function for fitting: y = m*x + b
from scipy.optimize import curve_fit
import numpy as np
from sklearn.metrics import r2_score

# Define linear function
def linear_func(x, a, b):
    return a * x + b

# Fit the data and get covariance matrix
popt_xx, pcov_xx = curve_fit(linear_func, strain, stress_xx)
popt_yy, pcov_yy = curve_fit(linear_func, strain, stress_yy)
popt_zz, pcov_zz = curve_fit(linear_func, strain, stress_zz)

# Extract standard errors (square root of diagonal elements of covariance matrix)
slope_xx_error = np.sqrt(pcov_xx[0, 0])
slope_yy_error = np.sqrt(pcov_yy[0, 0])
slope_zz_error = np.sqrt(pcov_zz[0, 0])

# Get fitted values
fitted_xx = linear_func(strain, *popt_xx)
fitted_yy = linear_func(strain, *popt_yy)
fitted_zz = linear_func(strain, *popt_zz)

residuals_xx = stress_xx - fitted_xx
residuals_zz = stress_zz - fitted_zz

# Compute R² values
r2_xx = r2_score(stress_xx, fitted_xx)
r2_yy = r2_score(stress_yy, fitted_yy)
r2_zz = r2_score(stress_zz, fitted_zz)

# Print results including errors
print(f"Fitted Line for Stress_xx: y = ({popt_xx[0]:.7f} ± {slope_xx_error:.7f}) * x + {popt_xx[1]:.5f}, R² = {r2_xx:.7f}")
print(f"Fitted Line for Stress_yy: y = ({popt_yy[0]:.7f} ± {slope_yy_error:.7f}) * x + {popt_yy[1]:.5f}, R² = {r2_yy:.7f}")
print(f"Fitted Line for Stress_zz: y = ({popt_zz[0]:.7f} ± {slope_zz_error:.7f}) * x + {popt_zz[1]:.5f}, R² = {r2_zz:.7f}")

# Create figure with two rows (top: stress plots, bottom: residuals), keeping them close together

fig, axes = plt.subplots(2, 2, figsize=(10, 6), gridspec_kw={'height_ratios': [3, 1.2]})

# Plot xx components
axes[0, 0].scatter(strain, stress_xx, label="Stress xx (data)", marker="o", s=1.5, alpha=0.7, color='blue')
axes[0, 0].plot(strain, fitted_xx, linestyle="-", color='black', linewidth = 0.7, label=f"Fit xx ")
axes[0, 0].set_ylabel("Stress (GPa)")
axes[0, 0].legend()


# Residual plot for xx (directly below main plot)
axes[1, 0].scatter(strain, residuals_xx, marker="o", s=1.5, color="blue", alpha=0.7)
axes[1, 0].axhline(0, linestyle="--", color="black", linewidth = 0.7)
axes[1, 0].set_xlabel("Strain")
axes[1, 0].set_ylabel("Residuals")

# Plot zz components
axes[0, 1].scatter(strain, stress_zz, label="Stress zz (data)", marker="o", s =1.5, alpha = 0.7, color = "blue")
axes[0, 1].plot(strain, fitted_zz, linestyle="-", color ='black', linewidth = 0.7,label=f"Fit zz ")
axes[0, 1].set_ylabel("Stress (GPa)")
axes[0, 1].legend()

# Residual plot for zz (directly below main plot)
axes[1, 1].scatter(strain, residuals_zz, marker="o", s=1.5, alpha = 0.7, color="blue")
axes[1, 1].axhline(0, linestyle="--", color="black",linewidth = 0.7)
axes[1, 1].set_xlabel("Strain")
axes[1, 1].set_ylabel("Residuals")

# Remove x-axis from main plot
axes[0, 1].set_xticklabels([])
axes[0, 1].set_xlabel("")
axes[0, 1].tick_params(axis='x', bottom=False)

# Remove x-axis from main plot
axes[0, 0].set_xticklabels([])
axes[0, 0].set_xlabel("")
axes[0, 0].tick_params(axis='x', bottom=False)

# Remove extra spacing and make plots tight
plt.subplots_adjust(hspace=0, wspace=0.25)


axes[0,0].set_xlim(-0.007,0.007)
axes[1,0].set_xlim(-0.007,0.007)
axes[0,1].set_xlim(-0.007,0.007)
axes[1,1].set_xlim(-0.007,0.007)

axes[0,0].set_ylim(-0.5,0.5)
axes[0,1].set_ylim(-1.3, 1.3)

axes[1,0].set_ylim(-0.02,0.03)
axes[1,1].set_ylim(-0.002,0.003)
axes[1,0].set_yticks([-0.01, 0, 0.01, 0.02])
axes[1,1].set_yticks([-0.002, 0, 0.002])

axes[1, 0].xaxis.set_major_locator(ticker.MaxNLocator(5))  # Residual plot for xx
axes[1, 1].xaxis.set_major_locator(ticker.MaxNLocator(5))  # Residual plot for zz


plt.show()
