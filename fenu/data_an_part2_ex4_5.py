#Ex 4-5 Data Analysis

import numpy as np 
import matplotlib.pylab as plt
import glob
import os

from cr_func import cr_spectrum, exp_particles

exposure = 18000 #km^2 sr yr #k_euso
telescope = "Auger"

#------------------------------------------------
#Ex 4.1
#Above 10^20 eV efficiency can be considered 100%

goal = 50 #cr

#Exp cr above 10^20 eV for k-euso (per year)
E_min = 20.0 # log scale [eV]
E_max = 22.0 # log scale [eV] - Default

part_p_yr = exp_particles(E_min, E_max, exposure, telescope) #exp cr above  10^20 eV per year

years = goal / part_p_yr

print(f"Needed years to collect 50 CRs above 10^20 eV: {years:.2f}")


#------------------------------------------------
#Ex 4.2
#Particle detected are proportional to exposure (keeping all other parameters invariate)
goal2 = 200 #cr
E_min = np.log10(5.0e19) #log scale [eV]

part_p_yr2 = exp_particles(E_min, E_max, exposure, telescope) #exp cr above  5.0*10^19 eV per year

#Particle detected are proportional to exposure (keeping all other parameters invariate)
cum_exp = goal2 / part_p_yr2 * exposure

#Exposure is proportional to exp. time
years2 = cum_exp / exposure

print(f"Needed cumulative exposure to collect 200 CRs above 5.0*10^19 eV: {cum_exp:.2f}")

print(f"Needed years to collect 200 CRs above 5.0*10^19 eV: {years2:.2f}")


#------------------------------------------------
#Ex 5.1 - 5.2 - 5.3

E_min = 10.**20.3 #log scale [eV]
E_max = 10.**20.4  #log scale [eV]
exp_yr = 5.0 #years

cumulat_exp = exposure * exp_yr
fc_level = 2.44 #Feldman Cousins confidence intervals

j_limit = fc_level / cum_exp / (E_max - E_min)

#flux limit is  inversely proportional to the exposure
j_limit_10exp = j_limit / 10.

j_limit_0_1_exp =  j_limit * 10.

print(f"90% CL limits to the flux for the K-EUSO detector after 5 years: {j_limit:.2e}")
print(f"90% CL limits to the flux for the K-EUSO detector after 5 years (exposure is 10 times higher): {j_limit_10exp:.2e}")
print(f"90% CL limits to the flux for the K-EUSO detector after 5 years (exposure is 10 times lower): {j_limit_0_1_exp:.2e}")