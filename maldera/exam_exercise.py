#Exam exercise for Data Analysis course (Part 3)

import numpy as np
import matplotlib.pyplot as plt


#----------------------------------------------------------------------------
# Detector parameters
area = 1.0  # m^2 (per plane)
distance = 0.4  # m (between planes)
logE_min, logE_max = 2., 6.  # log10(E) energy range
nbins = 50
logE_bin = np.linspace(logE_min, logE_max, nbins)
plane_width = np.sqrt(area)  # Assuming square planes: 1 m x 1 m


#----------------------------------------------------------------------------
#Functions

# Efficiency function
def efficiency(logE):
    return 1. - (1. / (np.exp((logE - 2.) / 0.2) + 1.))


# Monte Carlo simulation function
def simulate_effective_area(theta_deg, Nev=1000000):
    theta_rad = np.radians(theta_deg)  # Convert to radians
    xvett, yvett, zvett = np.sin(theta_rad), 0, np.cos(theta_rad)  # Directional cosines

    # Generate random impact points on the first plane (z = 0)
    rng = np.random.default_rng()
    x = -plane_width / 2 + plane_width * rng.random(Nev)
    y = -plane_width / 2 + plane_width * rng.random(Nev)

    # Project impact points onto the second plane (z = 0.4 m)
    xd = x + (distance / zvett) * xvett
    yd = y + (distance / zvett) * yvett

    # Generate energies with E^-1 spectrum using inverse transform sampling
    u1 = rng.random(Nev)
    energies = 10**(logE_min + u1 * (logE_max - logE_min))

    # Energy mask using the hit/miss method
    u5 = rng.random(Nev)
    e_mask = u5 < efficiency(np.log10(energies))

    # Geometric mask: Check if particles hit both planes
    g_mask = (
        (x > -plane_width / 2) & (x < plane_width / 2) &
        (y > -plane_width / 2) & (y < plane_width / 2) &
        (xd > -plane_width / 2) & (xd < plane_width / 2) &
        (yd > -plane_width / 2) & (yd < plane_width / 2)
    )

    # Combine energy and geometric masks
    total_mask = e_mask & g_mask

    #Bin data
    logE_sim = np.log10(energies)
    
    logE_sim_counts, _ = np.histogram(logE_sim, bins=logE_bin)
    logE_hits_counts, _ = np.histogram(logE_sim[total_mask], bins=logE_bin)

    # Calculate effective area
    projected_area = area * np.cos(theta_rad)  # projected area
    a_eff = projected_area * (logE_hits_counts / logE_sim_counts)

    return a_eff, logE_sim_counts, logE_hits_counts


#----------------------------------------------------------------------------

# Simulate for two angles
theta1, theta2 = 60., 25.
a_eff_60, logE_60, hits_60 = simulate_effective_area(theta1)
a_eff_25, logE_25, hits_25 = simulate_effective_area(theta2)


#----------------------------------------------------------------------------
# Plot results

#Efficiency plot
xE =  np.linspace(logE_min, logE_max, 1000)
plt.figure()
plt.plot(xE,efficiency(xE), label=r'$\epsilon(\log_{10}(E))$')
plt.xlabel('logE')
plt.ylabel('efficiency')
plt.legend()
plt.tight_layout()
plt.savefig("img/eff(E)_plot.pdf")


#Effective Area plot
bin_centers = 0.5 * (logE_bin[:-1] + logE_bin[1:])

fig1, ax1 = plt.subplots()
ax1.plot(bin_centers, a_eff_60, label=r'A$_{\text{eff}}$ - θ = 60°', color='blue')
ax1.set_xlabel(r'$\log_{10}(E)$')
ax1.set_ylabel(r'$A_{\text{eff}}(E)$ [m$^2$]')
ax1.legend()
fig1.tight_layout()
fig1.savefig("img/a_eff_60.pdf")

fig2, ax2 = plt.subplots()
ax2.plot(bin_centers, a_eff_25, label=r'A$_{\text{eff}}$ - θ = 25°', color='blue')
ax2.set_xlabel(r'$\log_{10}(E)$')
ax2.set_ylabel(r'$A_{\text{eff}}(E)$ [m$^2$]')
ax1.legend()
fig2.tight_layout()
fig2.savefig("img/a_eff_25.pdf")

#Simulated energy vs. hit energy
fig3, ax3 = plt.subplots()
ax3.plot(bin_centers, logE_60, label ="Simulated energy", color = "blue")
ax3.plot(bin_centers, hits_60, label = "Final energy", color = "orange")
ax3.set_xlabel(r'$\log_{10}(E)$')
ax3.set_ylabel("Counts")
ax3.set_title(r'$\theta = 60˚$')
ax3.legend()
fig3.tight_layout()
fig3.savefig("img/Energy_60.pdf")

fig4, ax4 = plt.subplots()
ax4.plot(bin_centers, logE_25, label ="Simulated energy", color = "blue")
ax4.plot(bin_centers, hits_25, label = "Final energy", color = "orange")
ax4.set_xlabel(r'$\log_{10}(E)$')
ax4.set_ylabel("Counts")
ax4.set_title(r'$\theta = 25˚$')
ax4.legend()
fig4.tight_layout()
fig4.savefig("img/Energy_25.pdf")