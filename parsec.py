# code to generate TAMS and base of the RGB files from Parsec models

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import glob

# models from https://people.sissa.it/~sbressan/CAF09_V1.2S_M36_LT/
# solar metallicity
files=glob.glob('Z0.017Y0.279/*')
# [M/H]~0.3
#files=glob.glob('Z0.03Y0.302/*')

plt.ion()
plt.clf()

f = open('tams_parsec.txt','w')
f2 = open('rgb_parsec.txt','w')

for i in range(0,50):
    data=ascii.read(files[i])
    print(i,files[i])
    
    #rad=np.log10(10**data['LOG_R']/6.96e10)
    rad=np.log10(np.sqrt(10**data['LOG_L'] * (10**data['LOG_TE']/5777.)**(-4.)))
    age=(data['AGE'])*1e-9
    
    um=np.where(data['PHASE'] > 4.)[0]
    plt.plot(data['LOG_TE'][um],rad[um],'.',color='black',ms=0.2)
    
    um=np.where((data['PHASE'] == 7.) & (age < 20.))[0]
    if (len(um) != 0):
        plt.plot(data['LOG_TE'][um],rad[um],'o',color='red')
        f.write('%12.5f %12.5f \n' % (10**data['LOG_TE'][um[0]],10**rad[um[0]]))
    
    um=np.where((data['PHASE'] == 8.) & (age < 20.))[0]
    if (len(um) != 0):
        plt.plot(data['LOG_TE'][um],rad[um],'o',color='green')
        f2.write('%12.5f %12.5f \n' % (10**data['LOG_TE'][um[0]],10**rad[um[0]]))
    
    #plt.draw()
    #input(':')
    
f.close()
f2.close()