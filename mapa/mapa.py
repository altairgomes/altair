from mpl_toolkits.basemap import Basemap
import numpy as np
import os
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, Angle

def coord_pack(coord):
    if type(coord) == SkyCoord:
        if not coord.isscalar:
            return coord
        return SkyCoord([coord.ra], [coord.dec], frame=coord.frame)
    a, b = [], []
    for i in coord:
        a.append(i.ra)
        b.append(i.dec)
    coord = SkyCoord(a,b, frame='fk5')
    return coord

def calcfaixa(vel, data, star, dist, ca, pa, tamanho, step, erro=None, ring=None, atm=None):
    vec = np.arange(0, int(8000/(np.absolute(vel.value))), step)
    g = np.sort(np.concatenate((vec,-vec[1:]), axis=0))
    latlon = {'clat':{'lon':[], 'lat':[], 'lab': [], 'x': [], 'y': [], 'labx': []}, 'lats': {'lon':[], 'lat':[], 'lon2':[], 'lat2':[], 'x': [], 'y': [], 'x2':[], 'y2':[]}}
    if not erro == None:
        latlon['erro'] = {'lon': [], 'lat': [], 'lon2':[], 'lat2':[]}
        err = erro*u.mas
        errd = (dist.to(u.km)*err.to(u.rad)).value*u.km
    if not ring == None:
        latlon['ring'] = {'lon': [], 'lat': [], 'lon2':[], 'lat2':[]}
    if not atm == None:
        latlon['atm'] = {'lon': [], 'lat': [], 'lon2':[], 'lat2':[]}
    pa = Angle(pa)
    pa.wrap_at('180d', inplace=True)
    if pa > 90*u.deg:
        paplus = pa - 180*u.deg
    elif pa < -90*u.deg:
        paplus = pa + 180*u.deg
    else:
        paplus = pa
    for delt in g:
        deltatime = delt*u.s
        datas1 = data + TimeDelta(deltatime)
        datas1.delta_ut1_utc = 0
        lon = star.ra - datas1.sidereal_time('mean', 'greenwich')
        m = Basemap(projection='ortho',lat_0=star.dec.value,lon_0=lon.value,resolution=None)
        a, b = m(lon.value, star.dec.value)
        a = a*u.m
        b = b*u.m
        dista = (dist.to(u.km)*ca.to(u.rad)).value*u.km
        ax = a + dista*np.sin(pa) + (deltatime*vel)*np.cos(paplus)
        by = b + dista*np.cos(pa) - (deltatime*vel)*np.sin(paplus)
        ax2 = ax - tamanho/2*np.sin(paplus)
        by2 = by - tamanho/2*np.cos(paplus)
        ax3 = ax + tamanho/2*np.sin(paplus)
        by3 = by + tamanho/2*np.cos(paplus)
        clon1, clat1 = m(ax.value, by.value, inverse=True)
        if delt == 0:
            latlon['clat']['cxy'] = [ax.value, by.value]
        if clon1 < 1e+30:
            latlon['clat']['lon'].append(clon1)
            latlon['clat']['lat'].append(clat1)
            latlon['clat']['lab'].append(datas1.iso)
        else:
            latlon['clat']['x'].append(ax.value)
            latlon['clat']['y'].append(by.value)
            latlon['clat']['labx'].append(datas1.iso)
        lon1, lat1 = m(ax2.value, by2.value, inverse=True)
        if lon1 < 1e+30:
            latlon['lats']['lon'].append(lon1) 
            latlon['lats']['lat'].append(lat1)
        else:
            latlon['lats']['x'].append(ax2.value) 
            latlon['lats']['y'].append(by2.value)
        lon2, lat2 = m(ax3.value, by3.value, inverse=True)
        if lon2 < 1e+30:
            latlon['lats']['lon2'].append(lon2) 
            latlon['lats']['lat2'].append(lat2)
        else:
            latlon['lats']['x2'].append(ax3.value) 
            latlon['lats']['y2'].append(by3.value)
        if not erro == None:
            ax2 = ax - errd*np.sin(paplus)
            by2 = by - errd*np.cos(paplus)
            ax3 = ax + errd*np.sin(paplus)
            by3 = by + errd*np.cos(paplus)
            lon1, lat1 = m(ax2.value, by2.value, inverse=True)
            if lon1 < 1e+30:
                latlon['erro']['lon'].append(lon1) 
                latlon['erro']['lat'].append(lat1)
            lon2, lat2 = m(ax3.value, by3.value, inverse=True)
            if lon2 < 1e+30:
                latlon['erro']['lon2'].append(lon2) 
                latlon['erro']['lat2'].append(lat2)
        if not ring == None:
            rng = ring*u.km
            ax2 = ax - rng*np.sin(paplus)
            by2 = by - rng*np.cos(paplus)
            ax3 = ax + rng*np.sin(paplus)
            by3 = by + rng*np.cos(paplus)
            lon1, lat1 = m(ax2.value, by2.value, inverse=True)
            if lon1 < 1e+30:
                latlon['ring']['lon'].append(lon1) 
                latlon['ring']['lat'].append(lat1)
            lon2, lat2 = m(ax3.value, by3.value, inverse=True)
            if lon2 < 1e+30:
                latlon['ring']['lon2'].append(lon2) 
                latlon['ring']['lat2'].append(lat2)
        if not atm == None:
            atmo = atm*u.km
            ax2 = ax - atmo*np.sin(paplus)
            by2 = by - atmo*np.cos(paplus)
            ax3 = ax + atmo*np.sin(paplus)
            by3 = by + atmo*np.cos(paplus)
            lon1, lat1 = m(ax2.value, by2.value, inverse=True)
            if lon1 < 1e+30:
                latlon['atm']['lon'].append(lon1) 
                latlon['atm']['lat'].append(lat1)
            lon2, lat2 = m(ax3.value, by3.value, inverse=True)
            if lon2 < 1e+30:
                latlon['atm']['lon2'].append(lon2) 
                latlon['atm']['lat2'].append(lat2)
    return latlon

def offset(datas, ca, pa, dist, vel, ob_off_ra, ob_off_de, st_off_ra, st_off_de):
    off_ra = ob_off_ra - st_off_ra
    off_de = ob_off_de - st_off_de
    dca = off_ra*np.sin(pa) + off_de*np.cos(pa)
    dt = -int(((off_ra*np.cos(pa) - off_de*np.sin(pa)).to(u.rad)*dist.to(u.km)/np.absolute(vel)).value)*u.s
    ca = ca + dca
    datas = datas + dt
    return ca, datas
