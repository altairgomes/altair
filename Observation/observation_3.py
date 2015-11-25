# coding=UTF-8

import numpy as np
import itertools
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5, EarthLocation, Angle, AltAz
from astropy.table import Table, Column, vstack
import os

######################################################################

ra_formatter = lambda x: "%07.4f" %x
dec_formatter = lambda x: "%06.3f" %x
alt_formatter = lambda x: "%3.1f" %x
dist_formatter = lambda x: "%3.0f" %x
int_formatter = lambda x: "%02d" %x
int2_formatter = lambda x: "%+03d" %x
sentinel = object()

def read(datafile, name_col=None, coord_col=None, comment_col=None, time_col=None, time_fmt='jd', skiprows=0):
    """
    Parameters
    ----------
    datafile : str
       Name of the file to read
    name_col : sequence, or list of numbers; optional
       Columns that will be stacked into the name of objects
    coord_col :  sequence, or list of numbers; optional
       Columns refered to the coordinates of the objects
    comment_col :  sequence, or list of numbers; optional
       Columns that will be stacked into the name of objects
    time_col :  sequence, or list of numbers; optional
       Columns that will be stacked into the name of objects
    time_fmt : str, optional
       Subformat for inputting string times
    """
    retornar = []
    if name_col:
        if type(name_col) not in [list, tuple, np.ndarray]:
            name_col = [name_col]
        nomes = np.loadtxt(datafile, skiprows=skiprows, usecols=(name_col), unpack=True, dtype ='S30', ndmin=1)
        if len(name_col) > 1:
            nome = nomes[0]
            for i in np.arange(len(name_col))[1:]:
                nome = np.core.defchararray.add(nome, ' ')
                nome = np.core.defchararray.add(nome, nomes[i])
        else:
            nome = nomes
        retornar.append(nome)
    if coord_col:
        if type(coord_col) not in [list, tuple, np.ndarray]:
            coord_col = [coord_col]
        coords = np.loadtxt(datafile, skiprows=skiprows, usecols=(coord_col), unpack=True, dtype ='S20', ndmin=1)
        coor = coords[0]
        for i in np.arange(len(coord_col))[1:]:
            coor = np.core.defchararray.add(coor, ' ')
            coor = np.core.defchararray.add(coor, coords[i])
        coord = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))
        retornar.append(coord)
    if comment_col:
        if type(comment_col) not in [list, tuple, np.ndarray]:
            comment_col = [comment_col]
        comments = np.loadtxt(datafile, skiprows=skiprows, usecols=(comment_col), unpack=True, dtype ='S30', ndmin=1)
        if len(comment_col) > 1:
            comment = comments[0]
            for i in np.arange(len(comment_col))[1:]:
                comment = np.core.defchararray.add(comment, ' ')
                comment = np.core.defchararray.add(comment, comments[i])
        else:
            comment = comments
        retornar.append(comment)
    if time_col:
        fmt_time = {'iso': ['S20'], 'jd': ['f8']}
        if type(time_col) not in [list, tuple, np.ndarray]:
            time_col = [time_col]
        times = np.loadtxt(datafile, skiprows=skiprows, usecols=(time_col), unpack=True, dtype =fmt_time[time_fmt][0], ndmin=1)
        if time_fmt == 'iso':  ####### iso nao deve funcionar por enquanto
            tim = times[0]
            a = len(time_col)
            len_iso = {2: [' '], 6: ['-', '-', ' ', ':',':']}
            for i in np.arange(len(time_col))[1:]:
                tim = np.core.defchararray.add(tim, len_iso[a][i-1]) 
                tim = np.core.defchararray.add(tim, coords[i])
            time = Time(tim, format=time_fmt, scale='utc')
        elif time_fmt == 'jd':
            time = Time(times, format='jd', scale='utc')
        retornar.append(time)
    return retornar
    
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
    
def close_obj(coord, size):
    coord = coord_pack(coord)
    ab, ba = np.meshgrid(np.arange(len(coord)), np.arange(len(coord)))
    sep = coord[ab].separation(coord[ba])
    c = np.where(sep < size)
    close = np.where(c[0] < c[1])
    pairs = np.vstack((c[0][close],c[1][close]))
    samefov = np.delete(np.arange(len(coord)), np.hstack((c[0][close],c[1][close])))
    samefov = samefov.reshape(len(samefov),1).tolist()
    n, m = np.unique(pairs[0], return_counts=True)
    y = np.in1d(pairs[0], n[np.where(m == 1)])
    n1, m1 = np.unique(pairs[1], return_counts=True)
    y1 = np.in1d(pairs[1], n1[np.where(m1 == 1)])
    samefov = samefov + pairs.T[y*y1].tolist()
    q = pairs.T[-(y*y1)]
    for z in np.unique(q.T[0]):
        b = q.T[1][np.where(q.T[0] == z)]
        combs = []
        for i in np.arange(len(b),0, -1):
            els = [[z] + list(x) for x in itertools.combinations(b, i)]
            combs.append(els)
        for i in combs:
            for a in i:
                d = False
                w = [list(x) for x in itertools.combinations(a,2)]
                if not np.all([i in q.tolist() for i in w]):
                    continue
                for v in samefov:
                    e = all(k in v for k in a)
                    d = bool(d + e)
                if d == False:
                    samefov = samefov + [a]
    return np.sort(samefov).tolist()

def midpoint_coord(coord, weighted=False, weight=[], ra_dec=False):
    """
    """
    i = coord_pack(coord)
    if weighted == False:
        x = (np.max(i.cartesian.x, axis=0) - np.min(i.cartesian.x, axis=0))/2 + np.min(i.cartesian.x, axis=0)
        y = (np.max(i.cartesian.y, axis=0) - np.min(i.cartesian.y, axis=0))/2 + np.min(i.cartesian.y, axis=0)
        z = (np.max(i.cartesian.z, axis=0) - np.min(i.cartesian.z, axis=0))/2 + np.min(i.cartesian.z, axis=0)
    else:
        if weight == []:
            weight = np.ones(i.shape)
        x = np.sum(i.cartesian.x*weight, axis=0)/np.sum(weight, axis=0)
        y = np.sum(i.cartesian.y*weight, axis=0)/np.sum(weight, axis=0)
        z = np.sum(i.cartesian.z*weight, axis=0)/np.sum(weight, axis=0)
    delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))
    alfa = np.arctan2(y, x)
    if ra_dec == True:
        return alfa,delta
    mean_coord = SkyCoord(alfa, delta, frame='fk5')
    dist_center = coord.separation(mean_coord)
    return mean_coord, dist_center

def precess(coord, timeut):
    """
    """
    while not timeut.isscalar:
        timeut=timeut[0]
    coord = coord_pack(coord)
    fk5_data = FK5(equinox=timeut)
    coord_prec = coord.transform_to(fk5_data)
    return coord_prec
    
def identify_min(coords, times, instants):
    """
    """
    a, b = np.meshgrid(np.arange(len(times)), np.arange(len(instants)))
    c = (times[a] - instants[b]).value
    d = np.absolute(c)
    e = np.argsort(d)
    f = e[:,0:2]
    ra, dec = np.array([]), np.array([])
    for j in np.arange(len(f)):
        if instants[j].iso in times[f[j]].iso:
            coor = coords[f[j][np.where(instants[j].iso == times[f[j]].iso)]][0]
        else:
            t = np.repeat([1./np.absolute((times[f[j]] - instants[j]).value)], len(coords[f[j]][0]), axis=0)
            coor = midpoint_coord(coords[f[j]], weighted=True, weight=t.transpose())[0]
        ra = np.concatenate((ra, coor.ra.to_string(unit=u.hourangle, precision=4)))
        dec = np.concatenate((dec, coor.dec.to_string(unit=u.deg, precision=3)))
    ra = ra.reshape((len(coords[f[j]][0]), len(ra)/len(coords[f[j]][0])))  ## teste1
    dec = dec.reshape((len(coords[f[j]][0]), len(dec)/len(coords[f[j]][0])))  ## teste1
    g = SkyCoord(ra.transpose(), dec.transpose(), unit=(u.hourangle, u.deg))
    return g
    
def mesh_coord(coord, time, ephem=None):
    """
    """
    ra = coord.ra
    if ephem:
        ra = np.append(ra, [0]*len(ephem.keys()))
    rs, ts = np.meshgrid(ra, time.sidereal_time('mean'))
    if ephem:
        ra = np.append(ra, [0]*len(ephem.keys()))
    return rs, ts

def sky_time(coord, time, rise_set=False, limalt=0*u.deg, site=EarthLocation(0.0, 0.0, 0.0), fuse=TimeDelta(0, format='sec', scale='tai')):
    """
    """
    if type(limalt) != u.quantity.Quantity:
        limalt = limalt*u.deg
    if time.isscalar == True:
        time = Time([time.iso], format='iso', scale='utc')
    coord = coord_pack(coord)
    timeut = time - fuse
    if len(time.shape) == 1:
        timeut = Time([[i] for i in timeut.jd], format='jd', scale='utc')
    timeut.delta_ut1_utc = 0
    timeut.location = site

#    modificar para caber ephemerides

    ra, ts = mesh_coord(coord, timeut)
    dif_h_sid = Angle(ra-ts)
    dif_h_sid.wrap_at('180d', inplace=True)
    dif_h_sol = dif_h_sid * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
    dif = TimeDelta(dif_h_sol.hour*u.h, scale='tai')
    culminacao = timeut + dif
    culminacao.delta_ut1_utc = 0
    culminacao.location = site
    if rise_set == True:
        hangle_lim = np.arccos((np.cos(90.0*u.deg-limalt) - np.sin(coord.dec)*np.sin(site.latitude)) / (np.cos(coord.dec)*np.cos(site.latitude)))
        tsg_lim = Angle(ra + hangle_lim)
        dtsg_lim = tsg_lim - culminacao.sidereal_time('mean')
        dtsg_lim.wrap_at(360 * u.deg, inplace=True)
        dtsg_lim_sol = dtsg_lim * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
        a = np.where(np.isnan(dtsg_lim_sol))
        dtsg_lim_sol[a] = Angle([48.0]*len(a[0])*u.hour)
        dtsg_np = TimeDelta((dtsg_lim_sol.hour*u.h))
        sunrise = culminacao - dtsg_np
        sunset = culminacao + dtsg_np
        if (site.latitude > 0*u.deg):
            alwaysup = np.where(coord.dec >= 90*u.deg - site.latitude + limalt)
            neverup = np.where(coord.dec <= -(90*u.deg - site.latitude - limalt))
        else:
            alwaysup = np.where(coord.dec <= -(90*u.deg + site.latitude + limalt))
            neverup = np.where(coord.dec >= 90*u.deg + site.latitude - limalt)
        culminacao = culminacao + fuse
        sunrise = sunrise + fuse
        sunset = sunset + fuse
        return culminacao, sunrise, sunset, alwaysup, neverup
    culminacao = culminacao + fuse
    return culminacao

def height_time(coord, time, time_left=False, limalt=0.0*u.deg, site=EarthLocation(0.0, 0.0, 0.0), fuse=TimeDelta(0, format='sec', scale='tai')):
    """
    """
    coord = coord_pack(coord)
    timeut = time - fuse
    timeut.location = site
    altaz = coord.transform_to(AltAz(obstime=timeut,location=site))
    if time_left == True:
        poente = sky_time(coord, time, rise_set=True, limalt=limalt, site=site, fuse=fuse)[2]
        time_rest = poente - time
        return altaz.alt, time_rest
    return altaz.alt

def instant_list(time_begin, time_end=None, time_step=TimeDelta(60*60, format='sec'), fmt='iso'):
    """
    """
    if not type(time_begin) == Time:
        time_begin = Time(time_begin, format=fmt, scale='utc')
    if not time_end:
        time_end = time_begin
    if not type(time_end) == Time:
        time_end = Time(time_end, format=fmt, scale='utc')
    if not type(time_step) == TimeDelta:
        time_step = TimeDelta(time_step, format='sec')
    c = np.arange(0,(time_end-time_begin).sec + time_step.sec,time_step.sec)
    tempo = time_begin + c*u.s
    return tempo
    
def resume_night(coord, time, name, limalt=0.0*u.deg, site=EarthLocation(0.0, 0.0, 0.0), fuse=TimeDelta(0, format='sec', scale='tai')):
    """
    """
    timeut = time - fuse
    coord = coord_pack(coord)
    coord_prec = precess(coord, timeut)
    culminacao, nascer, poente, alwaysup, neverup = sky_time(coord_prec, time, rise_set=True, limalt=limalt, site=site, fuse=fuse)
    ra, dec = coord.ra
    night = '\nRA: ' + ra + ', DEC: ' +  dec + ', Rise: ' + np.char.array(nascer.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Culmination: ' + \
np.char.array(culminacao.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Set: ' + np.char.array(poente.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + \
' TL, ' + np.char.array(name)
    return night

#################################################################################################################
    
class Observation(object):
    """
Location on the Earth for the observation instance.

Parameters
----------
fuse : number
    Fuse time of the site
latitude : sequence, or list of numbers; optional
    Columns that will be stacked into the name of objects
longitude :  sequence, or list of numbers; optional
    Columns refered to the coordinates of the objects
height :  sequence, or list of numbers; optional
    Columns that will be stacked into the name of objects
limheight :  sequence, or list of numbers; optional
    Columns that will be stacked into the name of objects
limdist : str, optional
    Subformat for inputting string times
    """

    def __init__(self, fuse = 0, latitude = 0.0, longitude = 0.0, height = 0.0, limheight=0.0, limdist=0.0):
        self.set_fuse(fuse)
        self.set_site(longitude, latitude, height)
        self.set_limheight(limheight)
        self.set_limdist(limdist)
        
    def set_site(self,longitude, latitude, height=0.0*u.m):
        """
Location on the Earth.

Initialization is first attempted assuming geocentric (x, y, z) coordinates
are given; if that fails, another attempt is made assuming geodetic
coordinates (longitude, latitude, height above a reference ellipsoid).
When using the geodetic forms, Longitudes are measured increasing to the
east, so west longitudes are negative. Internally, the coordinates are
stored as geocentric.

To ensure a specific type of coordinates is used, use the corresponding
class methods (`from_geocentric` and `from_geodetic`) or initialize the
arguments with names (``x``, ``y``, ``z`` for geocentric; ``lon``, ``lat``,
``height`` for geodetic).  See the class methods for details.


Notes
-----
This class fits into the coordinates transformation framework in that it
encodes a position on the `~astropy.coordinates.ITRS` frame.  To get a
proper `~astropy.coordinates.ITRS` object from this object, use the ``itrs``
property.
        """
        self.site = EarthLocation(longitude, latitude, height)

    def set_fuse(self, fuse):
        """
        """
        self.fuse = TimeDelta(fuse*3600, format='sec', scale='tai')
        
    def set_limheight(self, height):
        """
        """
        self.limheight = height*u.deg
        
    def set_limdist(self, size):
        """
        """
        self.limdist = size*u.arcmin

    def read(self, datafile, name_col=None, coord_col=None, comment_col=None, time_col=None, time_fmt='jd'):
        """
        """
        a = read(datafile, name_col, coord_col, comment_col, time_col, time_fmt)
        n = 0
        if name_col:
            self.names = a[n]
            n = n + 1
        else:
            self.names = ['']*np.shape(a)[1]
        if coord_col:
            self.coords = a[n]
            n = n + 1
        if comment_col:
            self.comments = a[n]
            n = n + 1
        else:
            self.comments = ['']*np.shape(a)[1]
        if time_col:
            self.times = a[n]
            n = n + 1
            
    def ephem(self, path='./ephemeris', coord_col=None, time_col=None, time_fmt='jd', skiprows=0):
        """
        """
        onlyeph = [ i for i in os.listdir(path) if i[-4:] == '.eph' ]
        ra, dec, body = np.array([]), np.array([]), [] ## teste1
        for i in onlyeph:
            a = read('{}/{}'.format(path,i), coord_col=coord_col, time_col=time_col, time_fmt=time_fmt, skiprows=skiprows)
            name = i[:-4]
            body.append(name) ## teste1
            ra = np.concatenate((ra, a[0].ra.to_string(unit=u.hourangle, precision=4))) ## teste1
            dec = np.concatenate((dec, a[0].dec.to_string(unit=u.deg, precision=3)))  ## teste1
        ra = ra.reshape((len(body), len(ra)/len(body)))  ## teste1
        dec = dec.reshape((len(body), len(dec)/len(body)))  ## teste1
        coord = SkyCoord(ra.transpose(), dec.transpose(), unit=(u.hourangle, u.deg))  ## teste1
        self.eph = {'name' : body, 'coord' : coord, 'time' : a[1]}  ## teste1
            
    def eph_moon(self, coord_col, time_col, ephem='Moon.eph', time_fmt='jd', skiprows=0):
        """
        """
        a = read(ephem, coord_col=coord_col, time_col=time_col, time_fmt=time_fmt, skiprows=skiprows)
        ra = a[0].ra.to_string(unit=u.hourangle, precision=4).reshape((1, len(a[0])))  ## teste1
        dec = a[0].dec.to_string(unit=u.deg, precision=3).reshape((1, len(a[0])))  ## teste1
        coord = SkyCoord(ra.transpose(), dec.transpose(), unit=(u.hourangle, u.deg))
        self.moon = {'coord' : coord, 'time' : a[1]}
        
    def __close_obj__(self):
        """
        """
        fov =  close_obj(self.coords, self.limdist)
        midcoord, midobj, midcomment = [], [], []
        for i in fov:
            if len(i) > 1:
                midcoord.append(midpoint_coord(self.coords[i])[0])
                a = ''
                for k in self.names[i][:-1]:
                    a = a + k + ' + '
                a = a + self.names[i][-1]
                midobj.append(a)
                b = ''
                for k in self.comments[i][:-1]:
                    b = b + k + ' + '
                b = b + self.comments[i][-1]
                midcomment.append(b)
            else:
                midcoord.append(self.coords[i[0]])
                midobj.append(self.names[i][0])
                midcomment.append(self.comments[i][0])
        midcoord = coord_pack(midcoord)
        self.samefov = {'coords': midcoord, 'names': np.array(midobj), 'comments': np.array(midcomment), 'fov': fov}
        
    def instant_list(self, time_begin, time_end=None, time_step=TimeDelta(60*60, format='sec')):
        """
        """
        self.instants = Time([[i] for i in instant_list(time_begin, time_end, time_step).jd], format='jd', scale='utc', location=self.site)
        self.instants.delta_ut1_utc = 0
        
    def plan(self, now=False, samefov=False):
        """
        """
        if hasattr(self, 'samefov'):
            coords = self.samefov['coords']
            names = self.samefov['names']
            comments = self.samefov['comments']
        elif hasattr(self, 'coords'):
            coords = self.coords
            names = self.names
            comments = self.comments
        if not hasattr(self, 'instants') or now == True:
            self.instant_list(Time.now() + self.fuse)
        instants = self.instants
        obs = Table(names=('objects', 'comments', 'RA_J2000_DEC', 'height', 'time_left', 'd_moon', 'culmination', 'time'),\
            dtype=(np.str, np.str, np.unicode, np.float, np.str, np.float, np.str, np.str))
        obs.meta = {'LT': instants[:,0], 'UT': (instants[:,0] - self.fuse)}
        if hasattr(self, 'moon'):
            moon_min = identify_min(self.moon['coord'], self.moon['time'], (instants - self.fuse)[:,0])
            coord_prec_moon = precess(moon_min, instants[0])
            alturam, time_restm = height_time(coord_prec_moon, instants, limalt=0*u.deg, time_left=True, site=self.site, fuse=self.fuse)
            obs.meta['moon_h'] = alturam.value.diagonal()
        if 'coords' in locals():
            coord_prec = precess(coords, instants[0])
            culmination, lixo, lixo2, alwaysup, neverup = sky_time(coord_prec, instants[0], limalt=self.limheight, rise_set=True, site=self.site, fuse=self.fuse)
            altura, time_rest = height_time(coord_prec, instants, limalt=self.limheight, time_left=True, site=self.site, fuse=self.fuse)
            if 'coord_prec_moon' in locals():
                ab, ba = np.meshgrid(np.arange(len(coords)), np.arange(len(coord_prec_moon)))
                distmoon = coords[ab].separation(coord_prec_moon[ba])
            for i in np.arange(len(instants)):
                t = Table()
                t['objects'] = Column(names)
                t['comments'] = Column(comments)
                t['RA_J2000_DEC'] = Column(coords.to_string('hmsdms', precision=4, sep=' '))
                t['height'] = Column(altura[i])
                t['time_left'] = Column(np.char.array([int_formatter(j) for j in time_rest[i].sec/3600.0]) + ':' + np.char.array([int_formatter(j) for j in (time_rest[i].sec - (time_rest[i].sec/3600.0).astype(int)*3600)/60]))
                if 'distmoon' in locals():
                    t['d_moon'] = Column(distmoon[i])
                t['culmination'] = Column(np.char.array(culmination[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0])
                t['time'] = instants[i,0].iso
                obs = vstack([obs, t[np.where(altura[i] > self.limheight)]])
        if hasattr(self, 'eph'):
            ephem_min = identify_min(self.eph[i]['coord'], self.eph[i]['time'], (instants - self.fuse)[:,0])
            coord_prec_eph = precess(ephem_min, instants[0])

##         continuar conferindo para ver se efemerides servem

            culminatione, lixo, lixo2, alwaysupe, neverupe = sky_time(coord_prec_eph, instants[0], limalt=self.limheight, rise_set=True, site=self.site, fuse=self.fuse)
            alturae, time_reste = height_time(coord_prec_eph, instants, limalt=self.limheight, time_left=True, site=self.site, fuse=self.fuse)
            alturae = alturae.diagonal();
            time_lefte = np.char.array([int_formatter(j) for j in time_reste.sec.diagonal()/3600.0]) + ':' + np.char.array([int_formatter(j) for j in (time_reste.sec.diagonal() - (time_reste.sec.diagonal()/3600.0).astype(int)*3600)/60])
            culmie = np.char.array(culminatione[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0]
            teph = Table([[i]*len(alturae), ephem_min.to_string('hmsdms', precision=4, sep=' '), alturae, time_lefte, culmie, instants[:,0].iso], names=('objects', 'RA_J2000_DEC', 'height', 'time_left', 'culmination', 'time'))                
            if 'coord_prec_moon' in locals():
                distmoon = coord_prec_eph.separation(coord_prec_moon)
                tephmoon = Column(distmoon, name='d_moon')
                teph.add_column(tephmoon, index=5)
            obs = vstack([obs, teph[np.where(alturae > self.limheight)]])
        obs['height'].format = '4.1f'
        obs['culmination'].format = '^10s'
        obs['culmination'].unit = 'LT'
        obs['d_moon'].format = '3.0f'
        obs['time_left'].format = '^8s'
        self.obs = obs

    
    def resume_night(self):
        """
        """
        self.night = resume_night(self.coords, self.instants[0], self.names, limalt=self.limheight, site=self.site, fuse=self.fuse)
            
    def create_plan(self, path='.', now=False):
        """
        """
        if not hasattr(self, 'instants') or now == True:
            self.instant_list(Time.now() + self.fuse)
        nome = '{}/Plano_{}.dat'.format(path, self.instants[0][0].iso.split(' ')[0])
        self.plan()
        #### imprime os dados de cada objeto para cada instante ####
        output = open(nome, 'w')
        output.write('Observational Plan to the night: {}\n\n'.format(self.instants[0][0].iso.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nMinimum height: {}\nField Size: {}\n\n'.format(self.site.latitude, self.site.longitude, self.limheight, self.limdist))
        output.write('Height: Height above the horizons (deg)\nTleft: Time left to reach minimum height (hh:mm)\n\n')
        for i in np.arange(len(self.instants)):
            obs = self.obs[np.where(self.obs['time'] == self.instants[i,0].iso)]
            obs.sort('time_left')
            output.write('\n---LT: {} (UT: {}) ----------------------------------------------------------------\n'.format(obs.meta['LT'][i].iso.rpartition(' ')[2].rpartition(':')[0], obs.meta['UT'][i].iso.rpartition(' ')[2].rpartition(':')[0]))
            if 'moon_h' in obs.meta:
                output.write('Moon height: {:4.1f}\n\n'.format(obs.meta['moon_h'][i]))
            for i in obs['RA_J2000_DEC', 'height', 'time_left', 'd_moon', 'culmination', 'objects', 'comments'].pformat(max_lines=-1, max_width=-1):
                output.write(i + '\n')
        self.resume_night()
        output.write('\n\n---Observability of the Targets----------------------------------------------------------------\n')
        for i in self.night[0]:
            output.write(i)
        output.close()
        
    def __region__(self, ra, dec, names, radius=10*u.arcsec):
        """
        """
        a = 'icrs\n'
        coords = SkyCoord(ra, dec, frame='fk5', unit=(u.hourangle, u.degree))
        for i in np.arange(len(coords)):
            a = a + 'circle({}, {}, {}) # text = '.format(coords[i].ra.deg, coords[i].dec.deg, (radius.to(u.deg)).value) + '{' + names[i] +'}\n'
        f = open('ds9.reg', 'w')
        f.write(a)
        f.close()
        
    def chart(self, size=None, server='eso', force_reg=False):
        """
        """
        dss = {'eso': '-dsseso', 'sao': '-dsssao'}
        if not size:
            size = self.limdist.to(u.arcmin).value
        if hasattr(self, 'obs'):
            p=0
        elif hasattr(self, 'samefov'):
            p =1
        else:
            p=2
        a = 0
        while a == 0:
            p = p + 1
            print p
            if p == 1:
                keys = self.instants[:,0].iso
                if len(keys) == 1:
                    m = 0
                else:
                    for l in np.arange(len(keys)):
                        print '{}: {} LT'.format(l, keys[l])
                    m = input('Choose the number of the date you want to show: ')
                key = keys[m]
                ra, dec = self.obs[key]['ra'], self.obs[key]['dec']
                names = np.concatenate((['List all objects'], self.obs[key]['names']))
            elif p == 2:
                if hasattr(self, 'samefov'):
                    ra, dec = self.samefov['coords'].to_string('hmsdms', precision=4, sep=' ')
                    names = np.concatenate((['List all objects'], self.samefov['names']))
                else:
                    continue
            elif p == 3:
                if hasattr(self, 'coords'):
                    ra, dec = self.coords.to_string('hmsdms', precision=4, sep=' ')
                    names = np.concatenate((['Exit'], self.names))
                else:
                    print "There is no coordinates to show"
                    return
            else:
                return
            print '\n'
            for i in np.arange(len(names)):
                print '{}: {}'.format(i, names[i])
            a = input('Choose the number of the target: ')
        if not os.path.isfile('ds9.reg') or force_reg == True:
            self.__region__(ra, dec, names[1:])
        os.system('ds9 {} size {} {} {} coord {} {} -region ds9.reg'.format(dss[server], size, size, dss[server], ra[a-1].replace(' ', ':'), dec[a-1].replace(' ', ':')))
        
    def show_text(self):
        """
        """
        k = self.instants[:,0].iso
        if len(k) == 1:
            i = 0
        else:
            for l in np.arange(len(k)):
                print '{}: {} LT'.format(l, k[l])
            i = input('Choose the number of the date you want to show: ')
        obs = self.obs[np.where(self.obs['time'] == k[i])]
        b = '\n---LT: {} (UT: {}) ----------------------------------------------------------------'.format(obs.meta['LT'][i].iso.rpartition(' ')[2].rpartition(':')[0], obs.meta['UT'][i].iso.rpartition(' ')[2].rpartition(':')[0])
        print b
        if 'moon_h' in obs.meta:
            print 'Moon height: {:+4.1f}'.format(obs.meta['moon_h'][i])
        obs.sort('time_left')
        print obs['RA_J2000_DEC', 'height', 'time_left', 'd_moon', 'culmination', 'objects', 'comments'].pprint(max_lines=-1, max_width=-1)
                
    def observe(self, edit=False):
        """
        """
        if edit == True:
            p=-1
        elif hasattr(self, 'obs'):
            p=0
        elif hasattr(self, 'samefov'):
            p =1
        else:
            p=2
        a = 0
        while a == 0:
            p = p + 1
            if p == 0:
                names = np.concatenate((self.observed))
            elif p == 1:
                key = self.obs.keys()[0]
                ra, dec = self.obs[key]['ra'], self.obs[key]['dec']
                names = np.concatenate((['List all objects'], self.obs[key]['names']))
            elif p == 2:
                if hasattr(self, 'samefov'):
                    coord = self.samefov['coords'].to_string('hmsdms', precision=4, sep=' ')
                    names = np.concatenate((['List all objects'], self.samefov['names']))
                else:
                    continue
            elif p == 3:
                if hasattr(self, 'coords'):
                    coord = self.coords.to_string('hmsdms', precision=4, sep=' ')
                    names = np.concatenate((['Exit'], self.names))
                else:
                    print "There is no coordinates to show"
                    return
            else:
                return
            print '\n'
            for i in np.arange(len(names)):
                print '{}: {}'.format(i, names[i])
            a = input('Digite o numero referente ao alvo: ')
        if edit == False:
            if not hasattr(self, 'observed'):
                self.observed = [names[a]]
            else:
                self.observed.append(names[a])
                
