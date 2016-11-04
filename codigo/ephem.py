import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from jplephem.spk import SPK
import astropy.constants as const
import os

#########################################################################################################

a = np.loadtxt('output.dat', unpack=True, dtype=np.float64)

t = a[0]
n = [801,802]
time = Time(t, format='jd', scale='tdb')

kernel = SPK.open('nep081.bsp')

for i in np.arange(len(n)):
    pos = a[3*i+1:3*i+4]
    vec = kernel[8,n[i]].compute(time.jd)[0:3] - kernel[8,899].compute(time.jd)[0:3]
    
    pos2 = (pos*u.AU).to(u.km).value
    d1 = np.sqrt(np.sum(pos2*pos2, axis=0))
    d2 = np.sqrt(np.sum(vec*vec, axis=0))
    
    
    l = ['X', 'Y', 'Z']
    for j in [0,1,2]:
        plt.plot(time.jyear, pos2[j]-vec[j], label=l[j])
    plt.plot(time.jyear, d1-d2, label='Distance')
    plt.xlabel('Time (years)')
    plt.ylabel(r'$\Delta$ Distance (km)')
    plt.title('(calculated - JPL) Object {}'.format(i+1))
    plt.legend(loc=3)
    plt.savefig('JPL_object_{}.png'.format(i+1), dpi=300)
    plt.clf()
