from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from observation_3 import read2

a = 'Nereida_eme.eph'
b = 'Nereida_jpl.eph'

c = read2(a, time_col=[0,1,2,3,4,5], time_fmt='iso', coord_col=[6,7], skiprows=7)
d = read2(b, time_col=[2], time_fmt='jd', coord_col=[4,5,6,7,8,9], skiprows=3)

eme_time = Time(c['times'], format='iso', scale='utc')
eme_coord = SkyCoord(c['coords'].transpose(), unit=(u.hourangle, u.deg))

jpl_time = Time(d['times'], format='jd', scale='utc')
jpl_coord = SkyCoord(d['coords'].transpose(), unit=(u.hourangle, u.deg))

dalfa = jpl_coord.ra*np.cos(jpl_coord.dec) - eme_coord.ra*np.cos(eme_coord.dec)
ddelta = jpl_coord.dec - eme_coord.dec

r = []
for i in np.arange(1992,2019,1):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

r = np.array(r)

plt.plot(jpl_time.jd - 2451544.5, dalfa.mas, label='Right Ascencion')
plt.plot(jpl_time.jd - 2451544.5, ddelta.mas, label='Declination')

plt.ylim(-500,500)
plt.title('Laurene - JUP300')
plt.xlabel('Time')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1992,2019,1)])
plt.ylabel('NEP081 - EMELYANOV (mas)')
plt.legend()
plt.axhline(0, color='black')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('JPL-EME_Nereida.png',dpi=100, bbox_inches='tight')
plt.clf()
