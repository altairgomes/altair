import numpy as np
import os
from multiprocessing import Pool

f = open('lista_naif', 'r')
obj = f.readlines()
f.close()

def gera_ephem(cod):
    cod = cod.strip()
    arq = 'teste_naif_{}'.format(cod)
    f = open(arq, 'w')
    f.write('874\n')
    f.write('{}\n'.format(cod))
    if cod[0] == '5':
        f.write('A\n')
    if cod == '999':
        f.write('J\n')
    f.write('de431\n')
    f.write('/home/altair/Documentos/teste/observation/jd_ephem.txt\n')
    f.close()
    os.system('./ephem_hv9l < {} > /home/altair/Documentos/teste/observation/Moon.eph'.format(arq, cod))
    os.system('rm {}'.format(arq))
    print 'feito: {}'.format(cod)
    
pool = Pool(processes=15)
pool.map(gera_ephem, obj)
