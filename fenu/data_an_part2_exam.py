#Exam Exercise Data Analysis - Part 2
#Lorenzo Cane - A.A 2024/2025
#assuming data in file are log_10(E) data where E is in eV

import numpy as np 
import matplotlib.pylab as plt
import glob
import os
from iminuit import Minuit
from iminuit.cost import LeastSquares


#---------------------------------------------------------------------------------------------
# Read the data file and define const
data = np.loadtxt('DataSet_Dec2024.txt')  
E_min = 2.0
E_max = 15.0
E_range = np.linspace(10.0**E_min, 10.0**E_max,10000)
step = 0.1
high_en_lim = 1.0e12 #eV

data_time = 6.3e7  #s
detector_surface = 2.5  #m^2
#detector_surface *= 1.0e-6 #km^2
diff_angle = 4. * np.pi

exposure = data_time * detector_surface * diff_angle #m^2 sr s 

energy_bins = np.arange(E_min, E_max, step)  # Logarithmic bins (ensure empty bin at the end)
#print(energy_bins)

#---------------------------------------------------------------------------------------------
# Bin the data (energy)
counts, _ = np.histogram(data, bins=energy_bins)
bin_centers = 10 ** ((energy_bins[:-1] + energy_bins[1:]) / 2)  # Center of each bin
bin_widths = 10 ** energy_bins[1:] - 10 ** energy_bins[:-1]

#---------------------------------------------------------------------------------------------
# Functions

def FC_limit(exposure, en_widht, N_lim = 2.44):
    #Calculate the flux limit acccording to Feldman - Cousins distribution
    # Exposure and energy bin widht must be provided
    
    j = N_lim / exposure / en_widht

    return j


def power_law(x, norm, index):
    return norm * x ** (-index)


#---------------------------------------------------------------------------------------------
# (a) Calculate the spectrum

spectrum_flux = np.divide(counts, (exposure)) / bin_widths

# Plot spectrum
plt.figure(figsize=(8, 6))
plt.scatter(bin_centers, spectrum_flux, label='Calculated Spectrum')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('E [eV]')
plt.ylabel(r'J(E) [$m^{-2} s^{-1} sr^{-1} eV^{-1} $] ')
plt.title('Spectrum')
#plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.savefig("img/Spectrum.pdf")

#---------------------------------------------------------------------------------------------
# (b) Expected particles in 10^9.3 to 10^9.4 eV bin

bin_index = np.where((energy_bins >= 9.3) & (energy_bins < 9.4))[0] #select correct energy bin
#print(bin_centers[bin_index])

bin_count_b = counts[bin_index[0]] if bin_index.size > 0 else 0
print(f"Particles in the bin (10^9.3 ; 10^9.4) eV: {bin_count_b}")

#---------------------------------------------------------------------------------------------
# (c) 90% FC upper limit for first and second empty bins

empty_bins_idx = np.where(counts == 0)[0] #select all empty bin indices
high_en_empty_bins = empty_bins_idx[bin_centers[empty_bins_idx] > high_en_lim] #select high en empty bins

upper_limits= []
for idx in high_en_empty_bins:
    em_bin_width = 10**(energy_bins[idx + 1]) - 10**(energy_bins[idx]) #empty en bins widht
    print(f"Empty bin: Energy range ({10**energy_bins[idx]:.1e} - {10**energy_bins[idx + 1]:.1e} ) eV, Bin width: {em_bin_width:.1e} eV")
    upper_limits.append(FC_limit(exposure, em_bin_width)) #upper limit

print(f"90% FC upper limits (1st empty bin): {upper_limits[0]: .4e}" + r' m^{-2} s^{-1} sr^{-1} eV^{-1}  ')
print(f"90% FC upper limits (2nd empty bin): {upper_limits[1]: .4e}" + r' m^{-2} s^{-1} sr^{-1} eV^{-1}  ')

#---------------------------------------------------------------------------------------------
# (d) Fit the data to a power law

error_flux = np.divide(np.sqrt(counts),exposure) / bin_widths # sqrt(N) as error on counts, exposure and energy bin w/o errors
error_flux[error_flux == 0] = 1.0e-30 #avoid 0 error


# Filter data for E (if ok)
mask_en = bin_centers > 10.0**8. 
mask_counts = counts != 0
mask = mask_en & mask_counts
E_filtered = bin_centers[mask]
flux_filtered = spectrum_flux[mask]
error_filtered = error_flux[mask]

assert len(bin_centers) == len(spectrum_flux)  #trying to avoid mismatch

l = LeastSquares(E_filtered, flux_filtered, error_filtered, model= power_law) #init least squares
m = Minuit(l, 1.0e-2,  2.5 , name=('A', "gamma"))
m.limits["A"] = (1.e-3, 1.e10)
m.limits["gamma"] = (0.5, 5.)
#m.errors["A"] = 1e-21
#m.errors["gamma"] = 0.1
m.migrad()
m.hesse()

print("Best fit parameters:\n")
for key, value, error in zip(m.parameters, m.values, m.errors):   #print fit results
    print(f"{key} : {value:.2e} +- {error:.2e}")

print(m.fmin)

# Plot data and fit
fitted_model = power_law(E_range, m.values["A"], m.values["gamma"])
plt.figure(figsize=(8, 6))
plt.scatter(bin_centers, spectrum_flux, label='Exp. Spectrum')
plt.plot(E_range, fitted_model, label='Fitted Power Law', color='red')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('E [eV]')
plt.ylabel(r'J(E) [$m^{-2} s^{-1} sr^{-1} eV^{-1} $] ')
plt.title('Power Law fitted Spectrum')
plt.legend()
plt.tight_layout()
#plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.savefig("img/Spectrum_filtered_fitted.pdf")

'''
#No filter
assert len(bin_centers) == len(spectrum_flux)  #trying to avoid mismatch

l = LeastSquares(bin_centers, spectrum_flux, error_flux, model= power_law) #init least squares
m = Minuit(l, 1.0e-2,  2.5 , name=('A', "gamma"))
m.limits["A"] = (1.e-10, 1.)
m.limits["gamma"] = (1., 5.)
#m.errors["A"] = 1e-21
#m.errors["gamma"] = 0.1
m.migrad()
m.hesse()

print("Best fit parameters:\n")
for key, value, error in zip(m.parameters, m.values, m.errors):   #print fit results
    print(f"{key} : {value:.2e} +- {error:.2e}")

print(m.fmin)

# Plot data and fit
fitted_model = power_law(E_range, m.values["A"], m.values["gamma"])
plt.figure(figsize=(8, 6))
plt.scatter(bin_centers, spectrum_flux, label='Exp. Spectrum')
plt.plot(E_range, fitted_model, label='Fitted Power Law', color='red')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('E [eV]')
plt.ylabel(r'J(E) [$m^{-2} s^{-1} sr^{-1} eV^{-1} $] ')
plt.title('Power Law fitted Spectrum')
plt.legend()
plt.tight_layout()
#plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.savefig("img/Spectrum_fitted.pdf")
'''

J_norm = m.values["A"]*(1.0e8)**(-m.values["gamma"])
print(f'Normalization at 10^8 eV : {J_norm:.2e}' + r' m^{-2} s^{-1} sr^{-1} eV^{-1}  ')
#---------------------------------------------------------------------------------------------
# (e) Years of acquisition for one particle in 10^14 to 10^14.1 eV

bin_width_e = 10 ** 14.1 - 10 ** 14.
bin_index_e = np.digitize(14, energy_bins)-1
if spectrum_flux[bin_index_e] > 0.0:
    flux_e = spectrum_flux[bin_index_e]
else:
    flux_e = FC_limit(exposure, bin_width_e)
time_required = 1. / (flux_e * detector_surface * bin_width_e) #time for detect one part in s
years_required = time_required / (365 * 24 * 3600) #s to years
print(f"Years required for one particle (in 10^14 to 10^14.1 eV): {years_required:.2f} yrs")

#---------------------------------------------------------------------------------------------
# (f) Impact of detector surface reduction by 30%

broken = 0.3 #broken part

new_surface = detector_surface * (1. - broken)
new_bin_count_b = counts[bin_index[0]] * new_surface / detector_surface
print(f"Particles in the bin with 30% reduced detector: {new_bin_count_b}")

#---------------------------------------------------------------------------------------------
# (g) Impact of 10% energy overestimation

#new_energy_bins = energy_bins - np.log10(1.1)  #dividing 1.1 or sub log_10(1.1)
counts_shifted, _ = np.histogram(data - np.log10(1.1), bins=energy_bins)
#counts_shifted = np.divide(counts_shifted, exposure) / bin_widths 
new_bin_count_g = counts_shifted[bin_index[0]] if bin_index.size > 0 else 0
print(f"Particles in the bin with 10% energy overestimation: {new_bin_count_g}")
