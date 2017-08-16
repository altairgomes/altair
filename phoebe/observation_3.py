# coding=UTF-8

import numpy as np
import itertools
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5, EarthLocation, Angle, AltAz
from astropy.table import Table, Column, vstack
import os

######################################################################

int_formatter = lambda x: "%02d" %x

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
            for i in np.arange(len(name_col), dtype=np.int8)[1:]:
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
        for i in np.arange(len(coord_col), dtype=np.int8)[1:]:
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
            for i in np.arange(len(comment_col), dtype=np.int8)[1:]:
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
            for i in np.arange(len(time_col), dtype=np.int8)[1:]:
                tim = np.core.defchararray.add(tim, len_iso[a][i-1]) 
                tim = np.core.defchararray.add(tim, coords[i])
            time = Time(tim, format=time_fmt, scale='utc')
        elif time_fmt == 'jd':
            time = Time(times, format='jd', scale='utc')
        retornar.append(time)
    return retornar
    
def read2(datafile, name_col=[], coord_col=[], comment_col=[], time_col=[], time_fmt='jd', skiprows=0):
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
    cols = name_col + coord_col + comment_col + time_col
    dados = np.loadtxt(datafile, skiprows=skiprows, usecols=cols, unpack=True, dtype ='S30', ndmin=2)
    retornar = {}
#        if type(name_col) not in [list, tuple, np.ndarray]:
#   coords         name_col = [name_col]
#        nomes = np.loadtxt(datafile, skiprows=skiprows, usecols=(name_col), unpack=True, dtype ='S30', ndmin=1)
    if len(name_col) > 0:
        nomes = np.arange(0,len(name_col), dtype=np.int8)
        nome = dados[0]
        for i in nomes[1:]:
            nome = np.core.defchararray.add(nome, ' ')
            nome = np.core.defchararray.add(nome, dados[i])
        retornar['names'] = nome
#    elif len(name_col) == 1:
#        retornar['names'] = dados[0]
    if len(coord_col) > 0:
#        if type(coord_col) not in [list, tuple, np.ndarray]:
#            coord_col = [coord_col]
        n = len(coord_col)
        coords = np.arange(0, n/2, dtype=np.int8) + len(name_col)
        ra = dados[coords[0]]
        dec = dados[coords[0]+n/2]
#        coor = dados[coords[0]]
        for i in coords[1:]:
            ra = np.core.defchararray.add(ra, ' ')
            ra = np.core.defchararray.add(ra, dados[i])
            dec = np.core.defchararray.add(dec, ' ')
            dec = np.core.defchararray.add(dec, dados[i+n/2])
#            coor = np.core.defchararray.add(coor, ' ')
#            coor = np.core.defchararray.add(coor, dados[i])
#        coord = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))
        retornar['coords'] = np.array([ra, dec])#coor
#    if len(comment_col) > 0:
#        if type(comment_col) not in [list, tuple, np.ndarray]:
#            comment_col = [comment_col]
    if len(comment_col) > 0:
        comments = np.arange(0,len(comment_col), dtype=np.int8) + len(name_col) + len(coord_col)
        comment = dados[comments[0]]
        for i in comments[1:]:
            comment = np.core.defchararray.add(comment, ' ')
            comment = np.core.defchararray.add(comment, comments[i])
        retornar['comments'] = comment
#    else:
#        comment = comments
#        retornar.append(comment)
    if len(time_col) > 0:
#        fmt_time = {'iso': ['S20'], 'jd': ['f8']}
        n = len(name_col) + len(coord_col) + len(comment_col)
        times = np.arange(0, len(time_col), dtype=np.int8) + n
#        if type(time_col) not in [list, tuple, np.ndarray]:
#            time_col = [time_col]
        if time_fmt == 'iso':  ####### iso nao deve funcionar por enquanto
            tim = dados[times[0]]
            a = len(time_col)
            len_iso = {2: [' '], 6: ['-', '-', ' ', ':',':']}
            for i in times[1:]:
                print a, i
                tim = np.core.defchararray.add(tim, len_iso[a][i-n-1]) 
                tim = np.core.defchararray.add(tim, dados[i])
#            time = Time(tim, format='iso', scale='utc')
        elif time_fmt == 'jd':
            tim = dados[times[0]].astype(np.float)
#            time = Time(dados[times[0]].astype(np.float), format='jd', scale='utc')
        retornar['times'] = tim
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
    ba, ab = np.indices((len(coord), len(coord)), dtype=np.int16)
    sep = coord[ab].separation(coord[ba])
    c = np.where(sep < size)
    close = np.where(c[0] < c[1])
    pairs = np.vstack((c[0][close],c[1][close]))
    samefov = np.delete(np.arange(len(coord), dtype=np.int16), np.hstack((c[0][close],c[1][close])))
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
        for i in np.arange(len(b),0, -1, dtype=np.int16):
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

def midpoint_coord(coord, center=False, weight=[], ra_dec=False, axis=0):
    """
    """
    i = coord_pack(coord)
    if center == True:
        x = (np.max(i.cartesian.x, axis=axis) + np.min(i.cartesian.x, axis=axis))/2
        y = (np.max(i.cartesian.y, axis=axis) + np.min(i.cartesian.y, axis=axis))/2
        z = (np.max(i.cartesian.z, axis=axis) + np.min(i.cartesian.z, axis=axis))/2
    else:
        if weight == []:
            weight = np.ones(i.shape, dtype=np.int8)
        x = np.sum(i.cartesian.x*weight, axis=axis)/np.sum(weight, axis=axis)
        y = np.sum(i.cartesian.y*weight, axis=axis)/np.sum(weight, axis=axis)
        z = np.sum(i.cartesian.z*weight, axis=axis)/np.sum(weight, axis=axis)
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
    b, a = np.indices((len(instants), len(times)), dtype=np.int16)
    c = (times[a] - instants[b]).value
    d = np.absolute(c)
    e = np.argsort(d)
    f = e[:,0:2]
    ra, dec = np.array([]), np.array([])
    for j in np.arange(len(f), dtype=np.int16):
        if instants[j].iso in times[f[j]].iso:
            coor = coords[f[j][np.where(instants[j].iso == times[f[j]].iso)]][0]
        else:
            t = np.ones((coords.shape[1], 2), dtype=np.int8)*1./np.absolute((times[f[j]] - instants[j]).value)
#            t2 = Time.now()
            coor = midpoint_coord(coords[f[j]], weight=t.transpose())[0]
#            t3 = Time.now()
#            print "midpoint: ", (t3-t2).sec
        ra = np.concatenate((ra, coor.ra.hourangle))
        dec = np.concatenate((dec, coor.dec.deg))
    ra = ra.reshape((len(ra)/len(coords[f[j]][0]), len(coords[f[j]][0])))  ## teste1
    dec = dec.reshape((len(dec)/len(coords[f[j]][0]), len(coords[f[j]][0])))  ## teste1
    g = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    return g
    
def mesh_coord(coord, time):
    """
    """
    ra = coord.ra
    time.delta_ut1_utc = 0
    if len(ra.shape) > 1:
        a, b = np.indices(ra.shape)
        ts = time.sidereal_time('mean')[a]
        return ra, ts
    rs, ts = np.meshgrid(ra, time.sidereal_time('mean'))
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
    ra, ts = mesh_coord(coord, timeut[:,0])
    dif_h_sid = Angle(ra-ts)
    dif_h_sid.wrap_at('180d', inplace=True)
    dif_h_sol = dif_h_sid * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
    dif = TimeDelta(dif_h_sol.hour*u.h, scale='tai')
    culminacao = timeut + dif
    culminacao.delta_ut1_utc = 0
    culminacao.location = site
    if (site.latitude > 0*u.deg):
        alwaysup = np.where(coord.dec >= 90*u.deg - site.latitude + limalt)
        neverup = np.where(coord.dec <= -(90*u.deg - site.latitude - limalt))
    else:
        alwaysup = np.where(coord.dec <= -(90*u.deg + site.latitude + limalt))
        neverup = np.where(coord.dec >= 90*u.deg + site.latitude - limalt)
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
        culminacao = culminacao + fuse
        sunrise = sunrise + fuse
        sunset = sunset + fuse
        return culminacao, sunrise, sunset, alwaysup, neverup
    culminacao = culminacao + fuse
    return culminacao, alwaysup, neverup

def height_time(coord, time, time_left=False, limalt=0.0*u.deg, site=EarthLocation(0.0, 0.0, 0.0), fuse=TimeDelta(0, format='sec', scale='tai')):
    """
    """
    coord = coord_pack(coord)
    timeut = time - fuse
    timeut.delta_ut1_utc = 0
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

class Coord(object):
    def __init__(self, name=[], ra=[], dec=[], datafile='', name_col=[], coord_col=[], mag_col=[], skiprows=0):
        if name.shape == ra.shape == dec.shape != []:
            self.name = name
            self.coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        else:
            a = read(datafile=datafile, name_col=name_col, coord_col=coord_col, comment_col=mag_col, skiprows=skiprows)
            n = 0
            if name_col:
                self.name = a[n]
                n = n + 1
            else:
                self.body = np.array(['']*np.shape(a)[0])
            if coord_col:
                self.coords = a[n]
                n = n + 1
            if mag_col:
                self.magnitudes = a[n]
                n = n + 1
            else:
                self.magnitudes = np.array(['']*np.shape(a)[0])
            
    def __read(self, datafile, name_col=[], coord_col=[], mag_col=[], time_fmt='jd', skiprows=0):
        """
        """
        a = read(datafile=datafile, name_col=name_col, coord_col=coord_col, comment_col=mag_col, time_fmt=time_fmt)
        n = 0
        if name_col:
            self.names = a[n]
            n = n + 1
        else:
            self.names = ['']*np.shape(a)[1]
        if coord_col:
            self.coords = a[n]
            n = n + 1
        if mag_col:
            self.magnitudes = a[n]
            n = n + 1
        else:
            self.comments = ['']*np.shape(a)[1]
            
    def midpoint_coord(self, weight=[], center=False, ra_dec=False, axis=0):
        """
        """
        coords = coord_pack(self.coords)
        if center == True:
            x = (np.max(coords.cartesian.x, axis=axis) + np.min(coords.cartesian.x, axis=axis))/2
            y = (np.max(coords.cartesian.y, axis=axis) + np.min(coords.cartesian.y, axis=axis))/2
            z = (np.max(coords.cartesian.z, axis=axis) + np.min(coords.cartesian.z, axis=axis))/2
        else:
            if weight == []:
                weight = np.ones(coords.shape, dtype=np.int8)
            x = np.sum(coords.cartesian.x*weight, axis=axis)/np.sum(weight, axis=axis)
            y = np.sum(coords.cartesian.y*weight, axis=axis)/np.sum(weight, axis=axis)
            z = np.sum(coords.cartesian.z*weight, axis=axis)/np.sum(weight, axis=axis)
        delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))
        alfa = np.arctan2(y, x)
        if ra_dec == True:
            return alfa,delta
        mean_coord = SkyCoord(alfa, delta, frame='fk5')
        dist_center = coords.separation(mean_coord)
        return mean_coord, dist_center
    
    def precess(self, timeut):
        """
        """
        while not timeut.isscalar:
            timeut=timeut[0]
        coord = coord_pack(self.coords)
        fk5_data = FK5(equinox=timeut)
        coord_prec = coord.transform_to(fk5_data)
        c = Coord(self.body, coord_prec.ra, coord_prec.dec)
        c.magnitudes = self.magnitudes
        return c
        
    def __mesh_coord(self, time):
        """
        """
        ra = self.coords.ra
        time.delta_ut1_utc = 0
        if len(ra.shape) > 1:
            a, b = np.indices(ra.shape)
            ts = time.sidereal_time('mean')[a]
            return ra, ts
        rs, ts = np.meshgrid(ra, time.sidereal_time('mean'))
        return rs, ts
    
    def sky_time(self, timeut, rise_set=False, limalt=0*u.deg, site=EarthLocation(0.0, 0.0, 0.0)):
        """
        """
        if type(limalt) != u.quantity.Quantity:
            limalt = limalt*u.deg
        if timeut.isscalar == True:
            timeut = Time([timeut.iso], format='iso', scale='utc')
        coord = coord_pack(self.coords)
        if len(timeut.shape) == 1:
            timeut = Time([[i] for i in timeut.jd], format='jd', scale='utc')
        timeut.delta_ut1_utc = 0
        timeut.location = site
        ra, ts = self.__mesh_coord(timeut[:,0])
        dif_h_sid = Angle(ra-ts)
        dif_h_sid.wrap_at('180d', inplace=True)
        dif_h_sol = dif_h_sid * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
        dif = TimeDelta(dif_h_sol.hour*u.h, scale='tai')
        culminacao = timeut + dif
        culminacao.delta_ut1_utc = 0
        culminacao.location = site
        if (site.latitude > 0*u.deg):
            alwaysup = np.where(coord.dec >= 90*u.deg - site.latitude + limalt)
            neverup = np.where(coord.dec <= -(90*u.deg - site.latitude - limalt))
        else:
            alwaysup = np.where(coord.dec <= -(90*u.deg + site.latitude + limalt))
            neverup = np.where(coord.dec >= 90*u.deg + site.latitude - limalt)
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
            culminacao = culminacao
            sunrise = sunrise
            sunset = sunset
            return culminacao, sunrise, sunset, alwaysup, neverup
        culminacao = culminacao
        return culminacao, alwaysup, neverup
    
    def height_time(self, timeut, time_left=False, limalt=0.0*u.deg, site=EarthLocation(0.0, 0.0, 0.0)):
        """
        """
        coord = coord_pack(self.coords)
        timeut.delta_ut1_utc = 0
        timeut.location = site
        altaz = coord.transform_to(AltAz(obstime=timeut,location=site))
        if time_left == True:
            poente = self.sky_time(timeut, rise_set=True, limalt=limalt, site=site)[2]
            time_rest = poente - timeut
            return altaz.alt, time_rest
        return altaz.alt
    
    def resume_night(self, timeut):
        """
        """
        coord = coord_pack(self.coords)
        coord_prec = self.precess(timeut)
        culminacao, nascer, poente, alwaysup, neverup = coord_prec.sky_time(timeut, rise_set=True, limalt=self.limheight, site=self.site)
        night = Table()
        night['Objects'] = Column(self.names)
        night['Comments'] = Column(self.magnitudes)
        night['RA_J2000_DEC'] = Column(coord.to_string('hmsdms', precision=4, sep=' '))
        night['Rise'] = Column(np.char.array(nascer[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0], format='^10s', unit='TL')
        night['Culmination'] = Column(np.char.array(culminacao[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0], format='^10s', unit='TL')
        night['Set'] = Column(np.char.array(poente[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0], format='^10s', unit='TL')
#        night = '\nRA: ' + ra + ', DEC: ' +  dec + ', Rise: ' + np.char.array(nascer.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Culmination: ' + \
#np.char.array(culminacao.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Set: ' + np.char.array(poente.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + \
#' TL, ' + np.char.array(name)
        return night

#################################################################################################################

class Ephemeris(object):
#__init__
#__read__
#identify_min

    def __init__(self, name=[], ra=np.array([]), dec=np.array([]), time=[], path='./ephemeris', coord_col=[], time_col=[], time_fmt='jd', skiprows=0):
        if ra.shape == dec.shape == (len(time), len(name)) != []:
            if type(name) == str:
                name = np.array(name)
#            ra = np.array(ra).reshape((len(ra), 1))
#            dec = np.array(dec).reshape((len(dec), 1))
        else:
            files = [ i for i in os.listdir(path) if i[-4:] == '.eph' ]
            body, ra, dec, time = self.__read(path, files, coord_col=coord_col, time_col=time_col, time_fmt=time_fmt, skiprows=skiprows)
        self.body = body
        self.coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        self.times = Time(time, format=time_fmt, scale='utc')
        
    def __read(self, path, files, coord_col=[], time_col=[], time_fmt='jd', skiprows=0):
        """
        """
        ra, dec, body = np.array([]), np.array([]), [] ## teste1
        for i in files:
            a = read2('{}/{}'.format(path,i), coord_col=coord_col, time_col=time_col, time_fmt=time_fmt, skiprows=skiprows)
            name = i[:-4]
            body.append(name) ## teste1
            ra = np.concatenate((ra, a['coords'][0])) ## teste1
            dec = np.concatenate((dec, a['coords'][1]))  ## teste1
        ra = ra.reshape((len(body), len(ra)/len(body)))  ## teste1
        dec = dec.reshape((len(body), len(dec)/len(body)))  ## teste1
        return np.array(body), ra.transpose(), dec.transpose(), a['times']  ## teste1

    def identify_min(self, instants):
        """
        """
        b, a = np.indices((len(instants), len(self.times)), dtype=np.int16)
        c = (self.times[a] - instants[b]).value
        d = np.absolute(c)
        e = np.argsort(d)
        f = e[:,0:2]
        ra, dec = np.array([]), np.array([])
        for j in np.arange(len(f), dtype=np.int16):
            if instants[j].iso in self.times[f[j]].iso:
                coor = self.coords[f[j][np.where(instants[j].iso == self.times[f[j]].iso)]][0]
            else:
                t = np.ones((self.coords.shape[1], 2), dtype=np.int8)*1./np.absolute((self.times[f[j]] - instants[j]).value)
                coor = midpoint_coord(self.coords[f[j]], weight=t.transpose())[0]
            ra = np.concatenate((ra, coor.ra.hourangle))
            dec = np.concatenate((dec, coor.dec.deg))
        ra = ra.reshape((len(ra)/len(self.coords[f[j]][0]), len(self.coords[f[j]][0])))  ## teste1
        dec = dec.reshape((len(dec)/len(self.coords[f[j]][0]), len(self.coords[f[j]][0])))  ## teste1
#        g = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        print ra.shape, dec.shape, len(self.body), len(instants)
        return Ephemeris(body=self.body, ra=ra, dec=dec, time=instants)

#################################################################################################################
    
class Observation(object):
#__init__ ok
#set_site ok
#set_fuse ok
#set_limheight ok
#set_limdist ok
#read ok
#ephem ok
#eph_moon ok
#__close_obj__
#instant_list ok
#plan ok
#resume_night - falta ephem
#create_plan ok
#__region__
#chart
#show_text ok
# observe
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

    def read(self, datafile, name_col=[], coord_col=[], comment_col=[], time_col=[], time_fmt='jd'):
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

    def read2(self, datafile, name_col=[], coord_col=[], comment_col=[], time_col=[], time_fmt='jd'):
        """
        """
        a = read2(datafile, name_col, coord_col, comment_col, time_col, time_fmt)
        if name_col:
            self.names = a['names']
        else:
            self.names = ['']*np.shape(a)[1]
        if coord_col:
            self.coords = a['coords']
        if comment_col:
            self.comments = a['comments']
        else:
            self.comments = ['']*np.shape(a)[1]
        if time_col:
            self.times = a['times']
            
    def ephem(self, path='./ephemeris', coord_col=None, time_col=None, time_fmt='jd', skiprows=0):
        """
        """
        onlyeph = [ i for i in os.listdir(path) if i[-4:] == '.eph' ]
        ra, dec, body = np.array([]), np.array([]), [] ## teste1
        for i in onlyeph:
            a = read2('{}/{}'.format(path,i), coord_col=coord_col, time_col=time_col, time_fmt=time_fmt, skiprows=skiprows)
            name = i[:-4]
            body.append(name) ## teste1
            ra = np.concatenate((ra, a['coords'][0])) ## teste1
            dec = np.concatenate((dec, a['coords'][1]))  ## teste1
        ra = ra.reshape((len(body), len(ra)/len(body)))  ## teste1
        dec = dec.reshape((len(body), len(dec)/len(body)))  ## teste1
        coord = SkyCoord(ra.transpose(), dec.transpose(), unit=(u.hourangle, u.deg))  ## teste1
        self.eph = {'name' : np.array(body), 'coord' : coord, 'time' : Time(a['times'], format=time_fmt, scale='utc')}  ## teste1
            
    def eph_moon(self, coord_col, time_col, ephem='Moon.eph', time_fmt='jd', skiprows=0):
        """
        """
        a = read2(ephem, coord_col=coord_col, time_col=time_col, time_fmt=time_fmt, skiprows=skiprows)
        ra = a['coords'][0].reshape((len(a['coords'][0]), 1))  ## teste1
        dec = a['coords'][1].reshape((len(a['coords'][1]), 1))  ## teste1
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        self.moon = {'coord' : coord, 'time' : Time(a['times'], format=time_fmt, scale='utc')}
        
    def __close_obj__(self):
        """
        """
        fov =  close_obj(self.coords, self.limdist)
        midcoord, midobj, midcomment = [], [], []
        for i in fov:
            if len(i) > 1:
                midcoord.append(midpoint_coord(self.coords[i], center=True)[0])
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
        obs = Table(names=('Objects', 'Comments', 'RA_J2000_DEC', 'Height', 'Time_left', 'D_moon', 'Culmination', 'time'),\
            dtype=(np.str, np.str, np.unicode, np.float, np.str, np.float, np.str, np.str))
        obs.meta = {'LT': instants[:,0], 'UT': (instants[:,0] - self.fuse)}
        if hasattr(self, 'moon'):
            moon_min = identify_min(self.moon['coord'], self.moon['time'], (instants - self.fuse)[:,0])
            coord_prec_moon = precess(moon_min, instants[0])
            alturam, time_restm = height_time(coord_prec_moon, instants, limalt=0*u.deg, time_left=True, site=self.site, fuse=self.fuse)
            obs.meta['moon_h'] = alturam[:,0].value
        if 'coords' in locals():
            coord_prec = precess(coords, instants[0])
            culmination, alwaysup, neverup = sky_time(coord_prec, instants, limalt=self.limheight, rise_set=False, site=self.site, fuse=self.fuse)
            altura, time_rest = height_time(coord_prec, instants, limalt=self.limheight, time_left=True, site=self.site, fuse=self.fuse)
            a, b = np.indices((altura.shape))
            t = Table()
            t['Objects'] = Column(names[b.reshape((b.size,))])
            t['Comments'] = Column(comments[b.reshape((b.size,))])
            t['RA_J2000_DEC'] = Column(coords[b.reshape((b.size,))].to_string('hmsdms', precision=4, sep=' '))
            t['Height'] = Column(altura.reshape((b.size,)))
            time_left = time_rest.sec.reshape((b.size,))
            t['Time_left'] = Column(np.char.array([int_formatter(j) for j in time_left/3600.0]) + ':' + np.char.array([int_formatter(j) for j in (time_left - (time_left/3600.0).astype(int)*3600)/60]))
            if 'coord_prec_moon' in locals():
                ba, ab = np.indices((len(coord_prec_moon), len(coords)), dtype=np.int16)
                distmoon = coords[ab].separation(coord_prec_moon)
                t['D_moon'] = Column(distmoon.reshape((b.size,)))
            culmi = culmination.iso.reshape((b.size,))
            t['Culmination'] = Column(np.char.array(culmi).rpartition(' ')[:,2].rpartition(':')[:,0])
            t['time'] = Column(instants[:,0].iso[a.reshape((a.size,))])
            obs = vstack([obs, t[np.where(altura.reshape((b.size,)) > self.limheight)]])
        t1 = Time.now()
        if hasattr(self, 'eph'):
            ephem_min = identify_min(self.eph['coord'], self.eph['time'], (instants - self.fuse)[:,0])
            t2 = Time.now()
            coord_prec_eph = precess(ephem_min, instants[0])
            culminatione, alwaysupe, neverupe = sky_time(coord_prec_eph, instants, limalt=self.limheight, rise_set=False, site=self.site, fuse=self.fuse)
            alturae, time_reste = height_time(coord_prec_eph, instants, limalt=self.limheight, time_left=True, site=self.site, fuse=self.fuse)
            t3 = Time.now()
            a, b = np.indices(ephem_min.shape)
            teph = Table()
            teph['Objects'] = Column(self.eph['name'][b.reshape((b.size,))])
#            t['comments'] = Column(comments)
            teph['RA_J2000_DEC'] = Column(ephem_min.to_string('hmsdms', precision=4, sep=' ').reshape((b.size,)))
            teph['Height'] = Column(alturae.reshape((b.size,)))
            time_lefte = time_reste.sec.reshape((b.size,))
            teph['Time_left'] = Column(np.char.array([int_formatter(j) for j in time_lefte/3600.0]) + ':' + np.char.array([int_formatter(j) for j in (time_lefte - (time_lefte/3600.0).astype(int)*3600)/60]))
            if 'coord_prec_moon' in locals():
                distmoon = coord_prec_eph.separation(coord_prec_moon)
                teph['D_moon'] = Column(distmoon.reshape((b.size,)))
            culmie = culminatione.iso.reshape((b.size,))
            teph['Culmination'] = Column(np.char.array(culmie).rpartition(' ')[:,2].rpartition(':')[:,0])
            teph['time'] = Column(instants[:,0].iso[a.reshape((a.size,))])
            obs = vstack([obs, teph[np.where(alturae.reshape((b.size,)) > self.limheight)]])
        t4 = Time.now()
#        print "min: ", (t2-t1).sec
#        print "culmi, alti: ", (t3-t2).sec
#        print "table", (t4-t3).sec
        obs['Height'].format = '4.1f'
        obs['Culmination'].format = '^10s'
        obs['Culmination'].unit = 'LT'
        obs['D_moon'].format = '3.0f'
        obs['Time_left'].format = '^8s'
        self.obs = obs

    def __resume_night__(self, time):
        """
        """
        timeut = time - self.fuse
        coord = coord_pack(self.coords)
        coord_prec = precess(coord, timeut)
        culminacao, nascer, poente, alwaysup, neverup = sky_time(coord_prec, time, rise_set=True, limalt=self.limheight, site=self.site, fuse=self.fuse)
        night = Table()
        night['Objects'] = Column(self.names)
        night['Comments'] = Column(self.comments)
        night['RA_J2000_DEC'] = Column(coord.to_string('hmsdms', precision=4, sep=' '))
        night['Rise'] = Column(np.char.array(nascer[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0], format='^10s', unit='TL')
        night['Culmination'] = Column(np.char.array(culminacao[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0], format='^10s', unit='TL')
        night['Set'] = Column(np.char.array(poente[0].iso).rpartition(' ')[:,2].rpartition(':')[:,0], format='^10s', unit='TL')
#        night = '\nRA: ' + ra + ', DEC: ' +  dec + ', Rise: ' + np.char.array(nascer.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Culmination: ' + \
#np.char.array(culminacao.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Set: ' + np.char.array(poente.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + \
#' TL, ' + np.char.array(name)
        return night
#        self.night = resume_night(self.coords, self.instants[0], self.names, limalt=self.limheight, site=self.site, fuse=self.fuse)
            
    def create_plan(self, path='.', now=False, sort='Time_left'):
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
        for i in np.arange(len(self.instants), dtype=np.int16):
            obs = self.obs[np.where(self.obs['time'] == self.instants[i,0].iso)]
            obs.sort(sort)
            output.write('\n---LT: {} (UT: {}) ----------------------------------------------------------------\n'.format(obs.meta['LT'][i].iso.rpartition(' ')[2].rpartition(':')[0], obs.meta['UT'][i].iso.rpartition(' ')[2].rpartition(':')[0]))
            if 'moon_h' in obs.meta:
                output.write('Moon height: {:4.1f}\n\n'.format(obs.meta['moon_h'][i]))
            for i in obs['RA_J2000_DEC', 'Height', 'Time_left', 'D_moon', 'Culmination', 'Objects', 'Comments'].pformat(max_lines=-1, max_width=-1):
                output.write(i + '\n')
        night = self.__resume_night__(self.instants[0])
        night.sort('Objects')
        output.write('\n\n' + '-'*100 + '\n')
        output.write('---Observability of the Targets----------------------------------------------------------------------\n')
        output.write('-'*100 + '\n')
        for i in night['RA_J2000_DEC', 'Rise', 'Culmination', 'Set', 'Objects','Comments'].pformat(max_lines=-1, max_width=-1):
            output.write(i + '\n')
        output.close()
        
    def __region__(self, ra, dec, names, radius=10*u.arcsec):
        """
        """
        a = 'icrs\n'
        coords = SkyCoord(ra, dec, frame='fk5', unit=(u.hourangle, u.degree))
        for i in np.arange(len(coords), dtype=np.int16):
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
                    for l in np.arange(len(keys), dtype=np.int16):
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
            for i in np.arange(len(names), dtype=np.int16):
                print '{}: {}'.format(i, names[i])
            a = input('Choose the number of the target: ')
        if not os.path.isfile('ds9.reg') or force_reg == True:
            self.__region__(ra, dec, names[1:])
        os.system('ds9 {} size {} {} {} coord {} {} -region ds9.reg'.format(dss[server], size, size, dss[server], ra[a-1].replace(' ', ':'), dec[a-1].replace(' ', ':')))
        
    def show_text(self, sort='Time_left'):
        """
        """
        k = self.instants[:,0].iso
        if len(k) == 1:
            i = 0
        else:
            for l in np.arange(len(k), dtype=np.int16):
                print '{}: {} LT'.format(l, k[l])
            i = input('Choose the number of the date you want to show: ')
        obs = self.obs[np.where(self.obs['time'] == k[i])]
        b = '\n---LT: {} (UT: {}) ----------------------------------------------------------------'.format(obs.meta['LT'][i].iso.rpartition(' ')[2].rpartition(':')[0], obs.meta['UT'][i].iso.rpartition(' ')[2].rpartition(':')[0])
        print b
        if 'moon_h' in obs.meta:
            print 'Moon height: {:+4.1f}'.format(obs.meta['moon_h'][i])
        obs.sort(sort)
        print obs['RA_J2000_DEC', 'Height', 'Time_left', 'D_moon', 'Culmination', 'Objects', 'Comments'].pprint(max_lines=-1, max_width=-1)
                
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
            for i in np.arange(len(names), dtype=np.int16):
                print '{}: {}'.format(i, names[i])
            a = input('Digite o numero referente ao alvo: ')
        if edit == False:
            if not hasattr(self, 'observed'):
                self.observed = [names[a]]
            else:
                self.observed.append(names[a])
                
