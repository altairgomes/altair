# coding=UTF-8

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
        lista_campos = []
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
                    for j in lista_campos:
                        e = all(k in j for k in a)
                        d = bool(d + e)
                    if d == True:
                        continue
                    mean_coord, dist_center = self.midpoint_coord(coord[a])
                    if np.isscalar(dist_center.value):
                        dist_center = [dist_center]
                    if all(k <= size/2 for k in dist_center):
                        lista_campos.append(a)
        return lista_campos
    
    def midpoint_coord(self, coord):
        """
        """
        i = self.coord_pack(coord)
        x = np.mean(np.cos(i.dec)*np.cos(i.ra))
        y = np.mean(np.cos(i.dec)*np.sin(i.ra))
        z = np.mean(np.sin(i.dec))
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
        
    def sky_time(self, coord, time, limalt=0*u.deg, site=None):
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
        hangle_lim = np.arccos((np.cos(90.0*u.deg-limalt) - np.sin(coord.dec)*np.sin(site.latitude)) / (np.cos(coord.dec)*np.cos(site.latitude)))
        tsg_lim = coord.ra + hangle_lim
        dtsg_lim = tsg_lim - culminacao.sidereal_time('mean')
        dtsg_lim.wrap_at(360 * u.deg, inplace=True)
        dtsg_lim_sol = dtsg_lim * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
        dtsg_np = TimeDelta(dtsg_lim_sol.hour*u.h)
        culminacao = culminacao + self.fuse
        sunrise = culminacao - dtsg_np
        sunset = culminacao + dtsg_np
        return culminacao, sunrise, sunset
        
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
            poente = self.sky_time(coord, time, limalt, site)[2]
            time_rest = poente - time
            return altura, time_rest
        return altura
        
    def plan(self, coord, time, obj=None, comment=None, limalt=0.0*u.arcmin, size=0.0*u.deg, site=None):
        if site == None:
            site = self.site
        if obj == None:
            obj = ['']*len(coord)
        if comment == None:
            comment = ['']*len(coord)
        comment1 = comment.copy()
        for i in np.arange(len(comment)):
            comment1[i] = '({})'.format(comment[i])
        if type(limalt) != u.quantity.Quantity:
            limalt = limalt*u.deg
        if type(size) != u.quantity.Quantity:
            size = size*u.arcmin
        timeut = time - self.fuse
        coord = self.coord_pack(coord)
        midcoord = np.copy(coord)
        midobj = obj.copy()
        midcomment = comment1.copy()
        lista_campos =  [np.array([0])]*len(midobj)
        if size > 0.0*u.arcmin:
            midcoord, midobj, midcomment = [], [], []
            lista_campos = self.close_obj(coord, size)
            for i in lista_campos:
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
        midcoord = self.coord_pack(midcoord)
        coord_prec = self.precess(midcoord, time)
        culmination = self.sky_time(coord_prec, time, site=site)[0]
        altura, time_rest = self.height_time(coord_prec, time, limalt=limalt, time_left=True)
        obs_tot = np.array([], dtype={'names': ('coord', 'height', 'time_left', 'obj', 'comment', 'culmination', 'tam_campo'),'formats': (object, 'f16', object, 'S100', 'S30', object, list)})
        for i in np.arange(len(altura)):
            a = [(midcoord[i], altura[i].value, time_rest[i], midobj[i][0], midcomment[i][0], culmination[i], lista_campos[i])]
            b = np.asarray(a, dtype={'names': ('coord', 'height', 'time_left', 'obj', 'comment', 'culmination', 'tam_campo'),'formats': (object, 'f16', object, 'S100', 'S30', object, list)})
            obs_tot = np.append(obs_tot,b)
        obs = np.sort(obs_tot[obs_tot['height']*u.deg >= limalt], order='time_left')
        a = '\n---LT: {} (UT: {})----------------------------------------------------------------\n'.format(time.iso.split(' ')[1][0:5], timeut.iso.split(' ')[1][0:5])
        for i in np.arange(len(obs)):
            a = a + '{} {}\n\tRA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n\tHeight: {:.1f}\n\tCulmination: {} LT\n\tTime left to reach min height: {:02d}:{:02d}\n'\
.format(obs['obj'][i], obs['comment'][i], obs['coord'][i].ra.hms.h, obs['coord'][i].ra.hms.m, obs['coord'][i].ra.hms.s, 
obs['coord'][i].dec.dms.d, np.absolute(obs['coord'][i].dec.dms.m), np.absolute(obs['coord'][i].dec.dms.s), obs['height'][i]*u.deg, obs['culmination'][i].iso.split(' ')[1][0:5],
int(obs['time_left'][i].sec/3600), int((obs['time_left'][i].sec - int(obs['time_left'][i].sec/3600)*3600)/60))
            if len(obs['tam_campo'][i]) > 1:
                for k in obs['tam_campo'][i]:
                    a = a + '\t\t{} {}\n\t\t  RA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n'.format(obj[k], comment1[k], coord[k].ra.hms.h,
coord[k].ra.hms.m, coord[k].ra.hms.s, coord[k].dec.dms.d, np.absolute(coord[k].dec.dms.m), np.absolute(coord[k].dec.dms.s))
            a = a + '\n'
        return a
    
    def resume_night(self, coord, time, name, site=None):
        if site == None:
            site = self.site
        timeut = time - self.fuse
        coord = self.coord_pack(coord)
        coord_prec = self.precess(coord, timeut)
        culminacao, nascer, poente = self.sky_time(coord_prec, timeut, limalt, site)
        a = '\n---Observability of the Targets----------------------------------------------------------------\n'
        for i in np.arange(len(nascer)):
            a = a + 'RA: {:02.0f} {:02.0f} {:07.4f}, DEC: {:+03.0f} {:02.0f} {:06.3f}, Sunrise: {} TL, Culmination: {} TL, Sunset: {} TL, {:10s}\n'\
.format(coord[i].ra.hms.h, coord[i].ra.hms.m, coord[i].ra.hms.s, coord[i].dec.dms.d, np.absolute(coord[i].dec.dms.m), np.absolute(coord[i].dec.dms.s),
nascer[i].iso.split(' ')[1][0:5], culminacao[i].iso.split(' ')[1][0:5], poente[i].iso.split(' ')[1][0:5], name[i])
        a = a + '\n'
        return a
    
    def create_plan(self, coord, obj, comment, horain, horafin, intinf, limalt=0.0*u.deg, size=0*u.arcmin, site=None, path='.'):
        """
        """
        if site == None:
            site = self.site
        tempoin = Time(horain, format='iso', scale='utc', location=site)
        tempofin = Time(horafin, format='iso', scale='utc', location=site) + TimeDelta(0.5, format='sec')
        intval = TimeDelta(intinf*60, format='sec')
        nome = '{}/Plano_{}'.format(path, horain.split(' ')[0])
        output = open(nome, 'w')
        output.write('Observational Plan to the night: {}\n\n'.format(horain.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nMinimun height: {}\nField Size: {}\n\n'.format(site.latitude, site.longitude, limalt, size))
        instante = tempoin
        while instante <= tempofin:
            #### imprime os dados de cada objeto para cada instante ####
            text = self.plan(coord, instante, obj, comment, site=site, limalt=limalt, size=size)
            output.write(text) 
            instante = instante + intval
        a = self.resume_night(coord, tempoin, obj)
        output.write(a)
        output.close()
        
#####################################################################

#arquivo = 'alvos.txt'				#### arquivo de alvos
#horain = '2015-03-08 18:00:00'			#### hora inicial local da observacao
#horafin = '2015-03-09 06:00:00'			#### hora final local da observacao
#intinf = 60					#### intervalo entre as informacoes (minutos)
#fuso = -3					#### fuso horario do local
#latitude = '-22 32 7.8'				#### latitude do local
#longitude = '314 25 2.5'			#### longitude do local
#altitude = 1864					#### altitude em metros
#limalt = 30.0					#### limite de altura para mostrar (graus)
#limdist = 11					#### limite de distancia para field-of-view (arcmin)

########################################################################

#observation = Observation(fuso, latitude, longitude, altitude)
#objs, coords, comments = observation.read(arquivo, [0,1,2], [3,4,5,6,7,8], [9])
#observation.create_plan(coords,objs, comments, horain, horafin, intinf, limalt=limalt, size=limdist)
#print observation.plan(coord, tempoin, objs, comments, limalt=30.0*u.deg, size=11.0*u.arcmin)

