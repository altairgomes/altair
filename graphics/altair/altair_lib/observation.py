# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
### BEGIN LICENSE
# Copyright (C) 2015 Altair Ramos altairgomesjr@gmail.com
# This program is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License version 3, as published 
# by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranties of 
# MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR 
# PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with this program.  If not, see <http://www.gnu.org/licenses/>.
### END LICENSE

import numpy as np
import itertools
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, Latitude, Longitude, FK5, EarthLocation, Angle

######################################################################

class Observation(object):
    """
    """

    def __init__(self, fuse = 0, latitude = 0.0, longitude = 0.0, height = 0.0):
        self.__set_fuse(fuse)
        self.__set_site(longitude, latitude, height)
        
    def __set_site(self,longitude, latitude, height=0.0*u.m):
        """
        """
        self.site = EarthLocation(longitude, latitude, height)

    def __set_fuse(self, fuse):
        """
        """
        self.fuse = TimeDelta(fuse*3600, format='sec', scale='tai')

    def read(self, datafile, name_col=None, coord_col=None, comment_col=None):
        """
        """
        nome, coord, comment, retornar = None, None, None, []
        if name_col != None:
            if type(name_col) != list and type(name_col) != tuple and type(name_col) != numpy.ndarray:
                name_col = [name_col]
            nomes = np.loadtxt(datafile, usecols=(name_col), unpack=True, dtype ='S30', ndmin=1)
            if len(name_col) > 1:
                nome = nomes[0]
                for i in np.arange(len(name_col))[1:]:
                    nome = np.core.defchararray.add(nome, ' ')
                    nome = np.core.defchararray.add(nome, nomes[i])
            else:
                nome = nomes
        if coord_col != None:
            if type(coord_col) != list and type(coord_col) != tuple and type(coord_col) != numpy.ndarray:
                coord_col = [coord_col]
            coords = np.loadtxt(datafile, usecols=(coord_col), unpack=True, dtype ='S20', ndmin=1)
            coor = coords[0]
            for i in np.arange(len(coord_col))[1:]:
                coor = np.core.defchararray.add(coor, ' ')
                coor = np.core.defchararray.add(coor, coords[i])
            coord = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))
        if comment_col != None:
            if type(comment_col) != list and type(comment_col) != tuple and type(comment_col) != numpy.ndarray:
                comment_col = [comment_col]
            comments = np.loadtxt(datafile, usecols=(comment_col), unpack=True, dtype ='S30', ndmin=1)
            if len(comment_col) > 1:
                comment = comments[0]
                for i in np.arange(len(comment_col))[1:]:
                    comment = np.core.defchararray.add(comment, ' ')
                    comment = np.core.defchararray.add(comment, comments[i])
            else:
                comment = comments
        return nome, coord, comment
        
    def coord_pack(self,coord):
        if type(coord) == SkyCoord:
            return coord
        a, b = [], []
        for i in coord:
            a.append(i.ra)
            b.append(i.dec)
        coord = SkyCoord(a,b, frame='fk5')
        return coord

    def close_obj(self, coord, size):
        """
        """
        coord = self.coord_pack(coord)
        if type(size) != u.quantity.Quantity:
            size = size*u.arcmin
        same_fov = []
        for idx, value in enumerate(coord):
            sep = coord[idx].separation(coord[idx+1:])
            close = [i + 1 + idx for i in np.arange(len(sep)) if sep[i] < size]
            combs = []
            for i in np.arange(len(close),0, -1):
                els = [[idx] + list(x) for x in itertools.combinations(close, i)]
                combs.append(els)
            combs.append([[idx]])
            for i in combs:
                for a in i:
                    d = False
                    for j in same_fov:
                        e = all(k in j for k in a)
                        d = bool(d + e)
                    if d == True:
                        continue
                    mean_coord, dist_center = self.midpoint_coord(coord[a])
                    if np.isscalar(dist_center.value):
                        dist_center = [dist_center]
                    if all(k <= size/2 for k in dist_center):
                        same_fov.append(a)
        return same_fov
    
    def midpoint_coord(self, coord):
        """
        """
        i = self.coord_pack(coord)
        x = np.mean(i.cartesian.x)
        y = np.mean(i.cartesian.y)
        z = np.mean(i.cartesian.z)
        delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))
        alfa = np.arctan2(y, x)
        mean_coord = SkyCoord(alfa, delta, frame='fk5')
        dist_center = coord.separation(mean_coord)
        return mean_coord, dist_center
        
    def precess(self, coord, time):
        """
        """
        timeut = time - self.fuse
        coord = self.coord_pack(coord)
        fk5_data = FK5(equinox=timeut)
        coord_prec = coord.transform_to(fk5_data)
        return coord_prec
        
    def sky_time(self, coord, time, rise_set=False, limalt=0*u.deg, site=None):
        """
        """
        if site == None:
            site = self.site
        if type(limalt) != u.quantity.Quantity:
            limalt = limalt*u.deg
        coord = self.coord_pack(coord)
        timeut = time - self.fuse
        timeut.delta_ut1_utc = 0
        timeut.location = site
        dif_h_sid = coord.ra - timeut.sidereal_time('mean')
        dif_h_sid.wrap_at('180d', inplace=True)
        dif_h_sol = dif_h_sid * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
        dif = TimeDelta(dif_h_sol.hour*u.h)
        culminacao = timeut + dif
        culminacao.delta_ut1_utc = 0
        if rise_set == True:
            hangle_lim = np.arccos((np.cos(90.0*u.deg-limalt) - np.sin(coord.dec)*np.sin(site.latitude)) / (np.cos(coord.dec)*np.cos(site.latitude)))
            tsg_lim = coord.ra + hangle_lim
            dtsg_lim = tsg_lim - culminacao.sidereal_time('mean')
            dtsg_lim.wrap_at(360 * u.deg, inplace=True)
            dtsg_lim_sol = dtsg_lim * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
            dtsg_np = TimeDelta(dtsg_lim_sol.hour*u.h)
            sunrise1 = culminacao - dtsg_np
            sunset1 = culminacao + dtsg_np
            sunrise, sunset = [], []
            for i in np.arange(len(sunset1)):
                if not np.isnan(hangle_lim[i]):
                    sunrise.append(sunrise1[i])
                    sunset.append(sunset1[i])
                    continue
                if (site.latitude > 0*u.deg and coord[i].dec >= 90*u.deg - site.latitude + limalt) or (site.latitude < 0*u.deg and coord[i].dec <= -(90*u.deg + site.latitude + limalt)):
                    sunrise.append('Always')
                    sunset.append('Always')
                else:
                    sunrise.append('Never')
                    sunset.append('Never')
            culminacao = culminacao + self.fuse
            return culminacao, sunrise, sunset
        culminacao = culminacao + self.fuse
        return culminacao
        
    def height_time(self, coord, time, time_left=False, limalt=0.0*u.deg, site=None):
        """
        """
        if site == None:
            site = self.site
        coord = self.coord_pack(coord)
        timeut = time - self.fuse
        timeut.location = site
        timeut.delta_ut1_utc = 0
        hourangle = timeut.sidereal_time('mean') - coord.ra
        distzen = np.arccos(np.sin(coord.dec)*np.sin(site.latitude) + np.cos(coord.dec)*np.cos(site.latitude)*np.cos(hourangle))
        altura = 90*u.deg - distzen
        if time_left == True:
            poente = self.sky_time(coord, time, rise_set=True, limalt=limalt, site=site)[2]
            time_rest = []
            for i in poente:
                if type(i) == Time:
                    time_rest.append(i - time)
                else:
                    time_rest.append(TimeDelta(2, format='jd'))
            return altura, time_rest
        return altura
        
    def __turn_time_array(self, time):
        if time.isscalar == False:
            return time
        return Time([time.iso])
        
    def time_list(self, time_begin, time_end, time_step):
        a = []
        tempo = time_begin
        while tempo <= time_end:
            a.append(tempo.iso)
            tempo = tempo + time_step
        return Time(a, format='iso', scale='utc')
        
    def plan(self, coord, time, obj, comment=None, limalt=0.0*u.deg, size=0.0*u.arcmin, site=None, same_fov=None):
        if site == None:
            site = self.site
        if obj == None:
            obj = ['']*len(coord)
        if comment == None:
            comment = ['']*len(coord)
        comment1 = comment
        for i in np.arange(len(comment)):
            comment1[i] = '({})'.format(comment[i])
        if type(limalt) != u.quantity.Quantity:
            limalt = limalt*u.deg
        if type(size) != u.quantity.Quantity:
            size = size*u.arcmin
        coord = self.coord_pack(coord)
        if size > 0.0*u.arcmin:
            midcoord, midobj, midcomment = [], [], []
            if same_fov == None:
                same_fov = self.close_obj(coord, size)
            for i in same_fov:
                if len(i) > 1:
                    midcoord.append(self.midpoint_coord(coord[i])[0])
                    a = ''
                    for k in obj[i][:-1]:
                        a  = a + k + ' + '
                    a = a + obj[i][-1]
                    midobj.append([a])
                    midcomment.append(np.array(['']))
                else:
                    midobj.append(obj[i])
                    midcoord.append(coord[i][0])
                    midcomment.append(comment1[i])
        else:
            midcoord = coord
            midobj = obj.copy()
            midcomment = comment1.copy()
            same_fov =  [np.array([0])]*len(midobj)
        midcoord = self.coord_pack(midcoord)
        time = self.__turn_time_array(time)
        coord_prec = self.precess(midcoord, time[0])
        culmination = self.sky_time(coord_prec, time[0], site=site)
        a = ''
        for tempo in time:
            altura, time_rest = self.height_time(coord_prec, tempo, limalt=limalt, time_left=True)
            obs_tot = np.array([], dtype={'names': ('coord', 'height', 'time_left', 'obj', 'comment', 'culmination', 'tam_campo'),'formats': (object, 'f16', 'S30', 'S100', 'S30', object, list)})
            for i in np.arange(len(altura)):
                if time_rest[i] < TimeDelta(1.5,format='jd'):
                    restante = '{:02d}:{:02d}'.format(int(time_rest[i].sec/3600), int((time_rest[i].sec - int(time_rest[i].sec/3600)*3600)/60))
                else:
                    restante = 'Always observed'
                arr = [(midcoord[i], altura[i].value, restante, midobj[i][0], midcomment[i][0], culmination[i], same_fov[i])]
                b = np.asarray(arr, dtype={'names': ('coord', 'height', 'time_left', 'obj', 'comment', 'culmination', 'tam_campo'),'formats': (object, 'f16', 'S30', 'S100', 'S30', object, list)})
                obs_tot = np.append(obs_tot,b)
            obs = np.sort(obs_tot[obs_tot['height']*u.deg >= limalt], order='time_left')
            a = a + '\n---LT: {} (UT: {}), N_objects={} ----------------------------------------------------------------\n'.format(tempo.iso.split(' ')[1][0:5], (tempo - self.fuse).iso.split(' ')[1][0:5], len(obs))
            for i in np.arange(len(obs)):
                a = a + '{} {}\n\tRA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n\tHeight: {:.1f} deg\n\tCulmination: {} LT\n\tTime left to reach min height: {}\n'\
.format(obs['obj'][i], obs['comment'][i], obs['coord'][i].ra.hms.h, obs['coord'][i].ra.hms.m, obs['coord'][i].ra.hms.s, 
obs['coord'][i].dec.dms.d, np.absolute(obs['coord'][i].dec.dms.m), np.absolute(obs['coord'][i].dec.dms.s), float(obs['height'][i]), obs['culmination'][i].iso.split(' ')[1][0:5],
obs['time_left'][i])
                if len(obs['tam_campo'][i]) > 1:
                    for k in obs['tam_campo'][i]:
                        a = a + '\t\t{} {}\n\t\t  RA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n'.format(obj[k], comment1[k], coord[k].ra.hms.h,
coord[k].ra.hms.m, coord[k].ra.hms.s, coord[k].dec.dms.d, np.absolute(coord[k].dec.dms.m), np.absolute(coord[k].dec.dms.s))
                a = a + '\n'
        return a
    
    def resume_night(self, coord, time, name, limalt=0*u.deg, site=None):
        if site == None:
            site = self.site
        timeut = time - self.fuse
        coord = self.coord_pack(coord)
        coord_prec = self.precess(coord, timeut)
        culminacao, nascer, poente = self.sky_time(coord_prec, time, rise_set=True, limalt=limalt, site=site)
        a = '\n---Observability of the Targets----------------------------------------------------------------\n'
        for i in np.arange(len(nascer)):
            if nascer[i] in ['Always', 'Never']:
                nas = nascer[i]
                poe = poente[i]
            else:
                nas = nascer[i].iso.split(' ')[1][0:5]
                poe = poente[i].iso.split(' ')[1][0:5]
            a = a + 'RA: {:02.0f} {:02.0f} {:07.4f}, DEC: {:+03.0f} {:02.0f} {:06.3f}, Rise: {} TL, Culmination: {} TL, Set: {} TL, {:10s}\n'\
.format(coord[i].ra.hms.h, coord[i].ra.hms.m, coord[i].ra.hms.s, coord[i].dec.dms.d, np.absolute(coord[i].dec.dms.m), np.absolute(coord[i].dec.dms.s),
nas, culminacao[i].iso.split(' ')[1][0:5], poe, name[i])
        a = a + '\n'
        return a
    
    def create_plan(self, coord, obj, comment, horain, horafin, intinf, limalt=0.0, size=0, site=None, path='.'):
        """
        """
        limalt = limalt*u.deg
        size = size*u.arcmin
        if site == None:
            site = self.site
        tempoin = Time(horain, format='iso', scale='utc', location=site)
        tempofin = Time(horafin, format='iso', scale='utc', location=site)
        intval = TimeDelta(intinf*60, format='sec')
        instante = self.time_list(tempoin,tempofin,intval)
        same_fov = self.close_obj(coord, size)
        nome = '{}/Plano_{}'.format(path, tempoin.iso.split(' ')[0])
        output = open(nome, 'w')
        output.write('Observational Plan to the night: {}\n\n'.format(tempoin.iso.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nMinimun height: {}\nField Size: {}\n\n'.format(site.latitude, site.longitude, limalt, size))
        #### imprime os dados de cada objeto para cada instante ####
        text = self.plan(coord, instante, obj, comment, site=site, limalt=limalt, size=size, same_fov=same_fov)
        output.write(text) 
        a = self.resume_night(coord, tempoin, obj, limalt=limalt)
        output.write(a)
        output.close()
        
#####################################################################

#arquivo = 'alvos'				#### arquivo de alvos
#horain = '2015-03-08 18:00:00'			#### hora inicial local da observacao
#horafin = '2015-03-09 06:00:00'			#### hora final local da observacao
#intinf = 60					#### intervalo entre as informacoes (minutos)
#fuso = -3					#### fuso horario do local
#latitude = '-22 32 7.8'				#### latitude do local
#longitude = '314 25 2.5'			#### longitude do local
#altitude = 1864					#### altitude em metros
#limalt = 0.0					#### limite de altura para mostrar (graus)
#limdist = 11					#### limite de distancia para field-of-view (arcmin)

#########################################################################

#observation = Observation(fuso, latitude, longitude, altitude)
#objs, coords, comments = observation.read(arquivo, [0], [1,2,3,4,5,6], [7])
#observation.create_plan(coords,objs, comments, horain, horafin, intinf, limalt=limalt, size=limdist)
#print observation.plan(coord, tempoin, objs, comments, limalt=30.0*u.deg, size=11.0*u.arcmin)

