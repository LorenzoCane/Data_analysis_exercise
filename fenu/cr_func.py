#Function used for cr exercise

import numpy as np


def cr_spectrum(E,telescope):

    j0 = 1.315e-18 #km^{-2} sr^{-1} yr^{-1} eV^{-1}

    if telescope == "Auger":
        p1 = 3.29
        p2 = 5.e18 #eV
        p3 = 2.51
        p4 = 1.3e19 #eV
        p5 = 3.05
        p6 = 4.6e19 #eV
        p7 = 5.1

    elif telescope == "TA":
        p1 = 3.23
        p2 = 5.4e18 #eV
        p3 = 2.63
        p4 = 1.8e19 #eV
        p5 = 2.92
        p6 = 7.1e19 #eV
        p7 = 5. 

    else : 
        raise ValueError("Telescope is not correct")
        return 0 

    E_val = 10**E
    term1 = j0 * (E_val / 10**18.5)**(-p1)
    term2 = (1 + (E_val / p2)**(1 / 0.05))**((p1 - p3) * 0.05)
    term3 = (1 + (E_val / p4)**(1 / 0.05))**((p3 - p5) * 0.05)
    term4 = (1 + (E_val / p6)**(1 / 0.05))**((p5 - p7) * 0.05)
    
    return term1 * term2 * term3 * term4






def exp_particles(Emin, Emax, exposure, telescope):
    """
    Calculate the expected number of events in the energy range [Emin, Emax], given the exposure of the instrument.
    
    Parameters:
        Emin (float): Minimum energy (log10 scale).
        Emax (float): Maximum energy (log10 scale).
        exposure (float): Exposure in km^2 sr yr.
    
    Returns:
        float: Calculated integral value.

    Exposure:
    mini_euso : 30000 km^2 sr yr
    k_euso : 18000 km^2 sr yr
    """
    # Perform the summation for the integral
    value = 0.0
    step = 0.1 # Step size in log10(E)
    for E in np.arange(Emin + 0.05, Emax, step):
        E_val_low = 10**(E - 0.05)
        E_val_high = 10**(E + 0.05)
        value += cr_spectrum(E,telescope) * (E_val_high - E_val_low) * exposure

    return value