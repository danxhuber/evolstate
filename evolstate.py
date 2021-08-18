import numpy as np
from astropy.io import ascii



def evolstate(teff, rad, logg, a=-0.68081514, b=2.3019458, c=1.5722924, d=16.269495, e=-37.498445,  
              path_tams='tams_parsec.txt', path_rgb='rgb_parsec.txt', path_mist='MIST_iso_5ae3ba3aad3cd.iso',
              teff_lim=5500.0, logg_lim=3.876, bin_teff=[5000.,5400.,6000.,6200.], bin_rad=[3.9,3.2,2.,1.],):
    
    """
    Simple function to assign an evolutionary state given the effective temperature, radius and/or surface
    gravity of a star (or many stars), where:
    0 == MS
    1 == SG
    2 == RGB
    3 == Binary (which can only be determined with teff & radius)
    
    Parameters
    ----------
    teff : np.array
        array of effective temperatures (K)
    rad : np.array
        array of stellar radii (in solar radii)
    logg : np.array
        array of surface gravities (dex)
        
    Returns
    -------
    cl : np.array
        array of stellar (evolutionary) classes i.e. 0-3
    
    """

    cl=np.zeros(len(teff))

    # teff+radius, see Berger et al. 2018
    if (rad[0] > -99):

        # Read in model grids
        sg_pc=ascii.read(path_tams)
        rgb_pc=ascii.read(path_rgb)
        mist=ascii.read(path_mist)
        
        # Keep only phase <= 3 (RGB)
        mask = np.where(mist['phase'] < 3)[0]
        mist_teff = np.array(10.0**mist['log_Teff'][mask])
        mist_rad = np.array(10.0**mist['log_R'][mask])
    
        # Empirical relation for dwarfs?
        mask = (mist_teff < 5029.0) & (mist_rad < 1.0)
        x = (mist_teff[mask]/4637.69 - 1.0)
        logl_lim =  a*x**0.0 + b*x**1.0 + c*x**2.0 + d*x**3.0 + e*x**4.0
        mist_rad[mask] = np.sqrt(1.55*(10.0**logl_lim) * (mist_teff[mask]/5777.)**(-4))
    
        # Define evol boundaries
        bin_cut = np.interp(rad, mist_rad, mist_teff)
        sg_cut = np.interp(teff, sg_pc['col1'], np.log10(sg_pc['col2']))
        rgb_cut = np.interp(teff, rgb_pc['col1'], np.log10(rgb_pc['col2']))    
    
        # Use boundaries to assign evolutionary states
        cl[(np.log10(rad) >= sg_cut) & (np.log10(rad) < rgb_cut)] = 1
        cl[(np.log10(rad) >= rgbcut)] = 2
        cl[(teff < bin_cut) & (teff > 3200.0)] = 3

    # teff+logg (no radius), see Huber et al. 2017
    if (logg[0] > -99):

        # Note to Dan: teff_lim variable is not actually used?
        glim_1 = np.interp(teff, bin_teff, bin_rad)
        glim = 1./4.671 * np.arctan((teff-6300.0)/-67.172) + 3.876
        
        cl[(logg < glim_1) & (logg < logg_lim)] = 2
        cl[(logg < glim) & (cl != 2)] = 1
        # For posterity, but not needed since the values start as 0
        cl[(logg >= glim) & (cl != 2)] = 0
    
    return cl
    


