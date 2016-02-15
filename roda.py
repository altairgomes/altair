import numpy as np
import os
from multiprocessing import Pool

################################

atual = os.getcwd()

f = open('lista_pastas')
pastas = f.readlines()
f.close()

def work(pasta):
    return

pool = Pool(processes=10)
pool.map(work, pastas)
###############################

a = open('PRAIA_astrometry_20_14.dat', 'w')
b = open('PRAIA_astrometry_20_11.dat', 'r')

c = b.readlines()

for i in c[0:6]:
    a.write(i)
a.write('output2                                           | extracted header data from fits images')
for i in c[6:42]:
    a.write(i)
for i in c[44:]:
    a.write(i)
    
    
from astropy.table import Table, vstack

f = open('lista_fits', 'r')
arqs = f.readlines()
f.close()

table = Table(names=('pasta', 'fits', 'stars'), dtype=(np.str, np.str, np.float))

for i in arqs:
    uc4 = i[:-5] + 'ucac4.red.xy'
    fil = i.split('/')
    if os.path.isfile(uc4):
        stars = np.loadtxt(uc4, usecols=(29,), ndmin=1)
        t = Table([[fil[0]], [fil[1]], [stars[0]]], names=('pasta', 'fits', 'stars'))
        table = vstack([table, t])
        continue
    t = Table([[fil[0]], [fil[1].strip()]], names=('pasta', 'fits'))
    table = vstack([table, t])
