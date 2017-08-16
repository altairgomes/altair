import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u

###############################################################

arquivo = 'search_results_ananke_OHP_final'
table = 'Ananke_ucac4_OHP.table_filtered'
output = 'Ananke_OHP_CDS.txt'
code = 'OH'
list_filt = {
'r cousins': 'R',
'v cousins': 'V'
}

###############################################################

f = open(arquivo, 'r')

dados = f.readlines()

f.close()

erros = np.loadtxt(table, usecols=(2, 3, 25), dtype={'names': ('ra', 'dec', 'time'), 'formats': ('f8', 'f8', 'f16')})
datasgeral = Time(erros['time'], format='jd', scale='utc')

dif = TimeDelta(60*60*6, format='sec')

g = open(output, 'w')

for i in dados:
    coord = SkyCoord(i[252:280], frame='icrs', unit=(u.hourangle,u.deg))
    time = np.float(i[303:320].strip())
    epoch = Time(time, scale='utc', format='jd')
    mag = np.float(i[79:86].strip())
    filtro = i[326:348].strip().lower()
    if filtro not in list_filt.keys():
        list_filt[filtro] = 'NAN'
        list_novo[filtro] = 'NAN'
    filt = list_filt[filtro]
    time_dif = np.absolute(datasgeral - epoch)
    idx = [i for i in np.arange(len(time_dif)) if time_dif[i] < dif]
    raerr = erros['ra'][idx]
    deerr = erros['dec'][idx]
    g.write(' {:02.0f} {:02.0f} {:07.4f} {:+03.0f} {:02.0f} {:06.3f} {:4.0f} {:4.0f} {:16.8f} {:4.1f} {:2s} {:2s}\n'\
.format(coord.ra.hms.h, coord.ra.hms.m, coord.ra.hms.s, coord.dec.dms.d, np.absolute(coord.dec.dms.m), np.absolute(coord.dec.dms.s), raerr[0], deerr[0], epoch.jd, mag, filt, code))

g.close()

print list_novo
