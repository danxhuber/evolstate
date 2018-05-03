import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# simple function to assign an evolutionary state given Teff and R/logg
# 0 = main sequence
# 1 = subgiant
# 2 = RGB
# 3 = main sequence binary (teff & R only)

def evolstate(teff,rad,logg):

    cl=np.zeros(len(teff))

    # Teff+Rad, see Berger et al. 2018
    if (rad[0] > -99):

        sg_pc=ascii.read('tams_parsec.txt')
        rgb_pc=ascii.read('rgb_parsec.txt')

        mist=ascii.read('MIST_iso_5ae3ba3aad3cd.iso')
        um=np.where(mist['phase'] < 3)[0]
        mistteff=10**mist['log_Teff'][um]
        mistrad=10**mist['log_R'][um]
    
        a=-0.68081514
        b=2.3019458
        c=1.5722924
        d=16.269495
        e=-37.498445
    
        um=np.where((mistteff < 5029.) & (mistrad < 1.))[0]
    
        x = (mistteff[um]/4637.69 - 1)
        logl_lim =  a + b*x + c*x**2 + d*x**3 + d*x**4
        mistrad[um]=np.sqrt(1.55*(10**logl_lim) * (mistteff[um]/5777.)**(-4))
    
        bincut=np.interp(rad,mistrad,mistteff)
    
        cl=np.zeros(len(teff))
        #cl[:]=1
    
        sgcut=np.interp(teff,sg_pc['col1'],np.log10(sg_pc['col2']))
        rgbcut=np.interp(teff,rgb_pc['col1'],np.log10(rgb_pc['col2']))    
    
        um=np.where((np.log10(rad) > sgcut) & (np.log10(rad) < rgbcut))[0]
        cl[um] = 1
    
        um=np.where((np.log10(rad) > rgbcut))[0]
        cl[um] = 2
    
        #um=np.where((np.log10(rad) > sgcut) & (np.log10(rad) < rgbcut) & (teff < 5300.) & (rad < 1.7))[0]
        #cl[um] = 3
    
        um=np.where((teff < bincut) & (teff > 3200.))[0]
        cl[um] = 3

    # Teff+logg, see Huber et al. 2017
    if (logg[0] > -99):
    
        tlim=5500.
        p=[6300.,-67.172,4.671,3.876]
        xl=[5000.,5400.,6000.,6200.]
        yl=[3.9,3.2,2.,1.]
        
        glim1=np.interp(teff,xl,yl)
        glim= 1./p[2] * np.arctan((teff-p[0])/p[1]) + p[3]
        
        u=np.where((logg < glim1) & (logg < p[3]))[0]
        cl[u]=2
        
        u=np.where((logg < glim) & (cl != 2))[0]
        cl[u]=1
        
        u=np.where((logg >= glim) & (cl != 2))[0]
        cl[u]=0
    
    return cl
    


