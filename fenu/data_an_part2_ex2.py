#Ex 2 Data Analysis

import numpy as np 
import matplotlib.pylab as plt
import glob
import os

from cr_func import cr_spectrum, exp_particles


#-----------------------------------------------------------------------
#Definitions and costants

Sim_l = 400.0 #km
FOV_l = 260.0 #km
nbins = 13    #number of bins

norm = (Sim_l / FOV_l)**2.0
#print(norm)
combined_file =  "combined.txt"

#telescope = "Auger"  #for flux calc
exposure = 18000 #km^2 sr yr k_euso



#-----------------------------------------------------------------------

#read events energy list file

energy_file = "OutputEnergy" #K-EUSO


en_sim_file = []
en_sim_num = []
energy_list = []
with open(energy_file, "r") as energies:
    for line in energies:
        en_sim_file.append(int(line.split()[0]))
        en_sim_num.append(int(line.split()[1]))
        en = (line.split()[2])
        energy_list.append(float(en)*1.0e6)
        #energy_list = np.multiply(energy_list, 1.0e6)

#print(energy_list)
#print(type(energy_list[0]))


#-----------------------------------------------------------------------
#plot hist of events energies

en_max = max(energy_list)
en_min = min(energy_list)
#print(en_max, en_min)
bins = np.logspace(np.log10(en_min), np.log10(en_max), nbins)

fig, ax = plt.subplots()
sim_count, sim_bins, sim_patches = ax.hist(energy_list, bins=bins)
ax.bar(bins[:-1], sim_count, width = np.diff(bins))
plt.xscale("log")
ax.set_xlabel("E [eV]")
ax.set_ylabel("Counts")

fig.savefig("img/energy_hist_k.pdf") #K-EUSO


#-----------------------------------------------------------------------
#read event trigger file

trigger_file = "OutputTrigger" #K-EUSO


tr_sim_file =[]
tr_sim_num = []
trigger_ev = []
with open(trigger_file, "r") as triggered:
    for line in triggered:
        tr_sim_file.append(int(line.split()[0]))
        tr_sim_num.append(int(line.split()[1]))
        trigger_ev.append(float(line.split()[2]))


#-----------------------------------------------------------------------
#select energy of events passing trigger process

trigger_en =[]

for sim, num, trig in zip(tr_sim_file, tr_sim_num, trigger_ev):
    for en_sim, en_num, ener in zip(en_sim_file, en_sim_num, energy_list):
        if en_sim == sim and en_num == num and trig!=0:
            trigger_en.append(ener)

#print(trigger_en)

#-----------------------------------------------------------------------
#calculate and plot efficiency


fig1, ax1= plt.subplots()
tr_count, tr_bins_edge = np.histogram(trigger_en, bins = bins)
#print(tr_count)

scaled_count = np.divide(tr_count, sim_count)
eff = np.multiply(scaled_count, norm)  #efficiency  
eff = np.clip(eff, 0.,1.)              #limit eff to [0,1] to avoid strange behaviour

threshold_index = np.argmax(eff >= 0.5)  # Identify the first bin where efficiency reaches or exceeds 0.5


ax1.bar(bins[:-1], eff, width=np.diff(bins), align='edge', color= "blue")


# Highlight the center of the bin with a vertical line
bin_center = (bins[threshold_index] + bins[threshold_index + 1]) / 2
ax1.axvline(bin_center, color="green", linestyle="--", label="50% Thrs.")
ax1.bar(bins[:-1], eff, width=np.diff(bins), align='edge')
ax1.axhline(0.5, color="r", linestyle="--", label="50% Eff.")


plt.xscale("log")
ax1.set_xlabel("E [eV]")
ax1.set_ylabel("Trigger Efficiency")
ax1.legend()
plt.tight_layout()
fig1.savefig("img/trigger_eff_k.pdf") #mini

print(f"50% Threshold level: {bin_center:.1e} eV")


#********************************************************************

parameters = [] #parameters for exp_particles function
for i in range(0, len(tr_bins_edge)-1):
    param = [np.log10(tr_bins_edge[i]), np.log10(tr_bins_edge[i+1]), exposure]
    parameters.append(param)


auger_expected_part = 0.0
ta_expected_part = 0.0
for p,eff in zip(parameters,eff):
    auger_coming_part =  exp_particles(p[0],p[1], p[2], "Auger")  #part count for each energy bins. Auger spectrum
    ta_coming_part =  exp_particles(p[0],p[1], p[2], "TA")  #part count for each energy bins. TA spectrum


    auger_expected_part += auger_coming_part * eff   #exp detected particles in each bin are summed. Auger spectrum
    ta_expected_part += ta_coming_part * eff   #exp detected particles in each bin are summed. TA spectrum


print(f"Number of events expected to trigger in 1 year (Auger spectrum) : {round(auger_expected_part):.1f}")
print(f"Number of events expected to trigger in 1 year (TA spectrum) : {round(ta_expected_part):.1f}")