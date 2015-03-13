# coding=UTF-8

import numpy as np
import itertools
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, Latitude, Longitude, FK5, EarthLocation, Angle

######################################################################

class Observation(object):

    def __init__(self):
        fuse = 0
        latitude = 0
        longitude = 0
        altitude = 0
        self.set_fuse(fuso)
        self.set_sitio(longitude, latitude, altitude)
        
    def read(self, arquivo):
        """
        """
        dados = np.loadtxt(arquivo, dtype={'names': ('objeto', 'data', 'hora', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'mag'), 
            'formats': ('S30', 'S30', 'S30', 'i4', 'i4', 'f8', 'i4', 'i4', 'f8', 'S30')})

        ### leitura data da ocultacao ###
        tempo = np.core.defchararray.add(dados['data'], [' '])
        tempo = np.core.defchararray.add(tempo, dados['hora'])
        time = Time(tempo, format='iso', scale='utc')
        
        ### leitura das coordenadas dos alvos ###
        alfa = np.core.defchararray.add(np.array(map(str, dados['afh'])), [' '])
        alfa = np.core.defchararray.add(alfa, np.array(map(str, dados['afm'])))
        alfa = np.core.defchararray.add(alfa, [' '])
        alfa = np.core.defchararray.add(alfa, np.array(map(str, dados['afs'])))
        delta = np.core.defchararray.add(np.array(map(str, dados['ded'])), [' '])
        delta = np.core.defchararray.add(delta, np.array(map(str, dados['dem'])))
        delta = np.core.defchararray.add(delta, [' '])
        delta = np.core.defchararray.add(delta, np.array(map(str, dados['afs'])))
        coor = np.core.defchararray.add(alfa, [' '])
        coor = np.core.defchararray.add(coor, delta)
        coordenadas = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))
        data = {'object': dados['objeto'], 'coord': coordenadas, 'time_occ': time, 'mag': dados['mag']}
        return data
        
    def set_sitio(self,longitude, latitude, altitude):
        """
        """
        self.sitio = EarthLocation(longitude, latitude, altitude)

    def set_fuse(self, fuso):
        """
        """
        self.utcoffset = TimeDelta(fuso*3600, format='sec', scale='tai')

    def close_obj(self, coord, limdist):
        """
        """
        tamcampo = Angle(limdist*u.arcmin)
        lista_campos = []
        for idx, value in enumerate(coord[:-1]):
            sep = coord[idx].separation(coord[idx+1:])
            close = [i + 1 + idx for i in np.arange(len(sep)) if sep[i] < tamcampo]
            combs = []
            for i in np.arange(len(close),0, -1):
                els = [list(x) for x in itertools.combinations(close, i)]
                combs.append(els)
            for i in combs:
                for l in i:
                    a = np.concatenate(([idx], l))
                    mean_coord, dist_center = self.midpoint_coord(coord[a])
                    if all(k <= tamcampo/2 for k in dist_center):
                        d = False
                        for j in lista_campos:
                            e = all(k in j for k in a)
                            d = bool(d + e)
                        if d == False:
                            lista_campos.append(a)
            d = False
            for i in lista_campos:
                e = idx in i
                d = bool(d + e)
            if d == False:
                lista_campos.append(np.array([idx]))
        return lista_campos
    
    def midpoint_coord(self, coords):
        """
        """
        x = np.mean([np.cos(i.dec)*np.cos(i.ra) for i in coords])
        y = np.mean([np.cos(i.dec)*np.sin(i.ra) for i in coords])
        z = np.mean([np.sin(i.dec) for i in coords])
        delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))*u.rad
        alfa = np.arctan(y/x)*u.rad
        mean_coord = SkyCoord(alfa, delta, frame='fk5')
        dist_center = coords.separation(mean_coord)
        return mean_coord, dist_center
        
    def precess(self, coord, time):
        """
        """
        fk5_data = FK5(equinox=time)
        coord_prec = coord.transform_to(fk5_data)
        return coord_prec
        
    def sky_time(self, coord, tempo, limalt=0*u.deg, sitio=self.sitio):
        """
        """
        tempo.delta_ut1_utc = 0
        tempo.location = sitio
        dif_h_sid = coord.ra - tempo.sidereal_time('mean')
        dif_h_sid.wrap_at('180d', inplace=True)
        dif_h_sol = dif_h_sid * (23 + 56/60 + 4.0916/3600) / 24
        dif = TimeDelta(dif_h_sol.hour*60*60, format='sec')
        culminacao = tempo + dif + self.utcoffset
        culminacao.delta_ut1_utc = 0
        hangle_lim = np.arccos((np.cos(90*u.deg-limalt) - np.sin(coord.dec)*np.sin(sitio.latitude)) / np.cos(coord.dec)*np.cos(sitio.latitude))
        tsg_lim = coord.ra + hangle_lim
        dtsg_lim = tsg_lim - culminacao.sidereal_time('mean')
        dtsg_lim.wrap_at(360 * u.deg, inplace=True)
        dtsg_lim_sol = dtsg_lim * (23 + 56/60 + 4.0916/3600) / 24
        dtsg_np = TimeDelta(dtsg_lim_sol.hour*60*60, format='sec')
        nascer = culminacao - dtsg_np
        poente = culminacao + dtsg_np
        return culminacao, nascer, poente
        
    def altura_time(self, coord, instante, sitio=self.sitio, limalt=0*u.deg):
        """
        """
        instante.location = sitio
        instante.delta_ut1_utc = 0
        hourangle = instante.sidereal_time('mean') - coord.ra
        distzen = np.arccos(np.sin(coord.dec)*np.sin(sitio.latitude) + np.cos(coord.dec)*np.cos(sitio.latitude)*np.cos(hourangle))
        altura = 90*u.deg - distzen
        poente = self.sky_time(coord, instante, limalt)[1]
        time_rest = poente - instante
        return altura, time_rest
        
    def create_plan(self, coord, obj, data, mag, horain, horafin, intinf, limalt=0*u.deg, limdist=0*u.deg sitio=self.sitio):
        """
        """
        tempoin = Time(horain, format='iso', scale='utc', location=sitio) - self.utcoffset
        tempofin = Time(horafin, format='iso', scale='utc', location=sitio) + TimeDelta(0.5, format='sec')  - self.utcoffset
        intval = TimeDelta(intinf*60, format='sec')
        coord = self.precess(coord, tempoin)
        culminacao, nascer, poente = self.sky_time(coord, tempoin, limalt, sitio)
        nome = 'Plano_{}'.format(horain.split(' ')[0])
        output = open(nome, 'w')
        output.write('Plano de observação da noite {}\n\n'.format(horain.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nLimite de altura: {} deg\nTamanho do campo: {} arcmin\n\n'.format(sitio.latitude, sitio.longitude, limalt, limdist))
        instante = tempoin
        while instante <= tempofin:
            altura, time_rest = self.altura_time(instante)
            #### imprime os dados de cada objeto para cada instante ####
            output.write('\n---LT: {} (UT: {})----------------------------------------------------------------\n'.format((instante + self.utcoffset).iso.split(' ')[1][0:5], instante.iso.split(' ')[1][0:5]))
            for idx, val in enumerate(altura):
                if altura[idx] > limalt*u.deg:
                    output.write('{} {} ({})\n\tRA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n\tAltura: {:.1f}\n\tCulminacao: {} LT\n\tTempo restante para limite: {:02d}:{:02d}\n\n'\
.format(obj[idx], data[idx], mag[idx], coord[idx].ra.hms.h, coord[idx].ra.hms.m, coord[idx].ra.hms.s,
coord[idx].dec.dms.d, coord[idx].dec.dms.m, coord[idx].dec.dms.s, altura[idx], culminacao[idx].iso.split(' ')[1][0:5],
int(time_rest[idx].sec/3600), int((time_rest[idx].sec - int(time_rest[idx].sec/3600)*3600)/60))) 
            instante = instante + intval
        output.write('\n---Observabilidade dos Alvos----------------------------------------------------------------\n')
        for idx, val in enumerate(nascer):
            output.write('{:10s} {}, RA: {:02.0f} {:02.0f} {:07.4f}, DEC: {:+03.0f} {:02.0f} {:06.3f}, Nascer: {} TL, Culminacao: {} TL, Poente: {} TL\n'\
.format(objeto[idx], data[idx], coord[idx].ra.hms.h, coord[idx].ra.hms.m, coord[idx].ra.hms.s,
coord[idx].dec.dms.d, np.absolute(coord[idx].dec.dms.m), np.absolute(coord[idx].dec.dms.s), nascer[idx].iso.split(' ')[1][0:5], culminacao[idx].iso.split(' ')[1][0:5],
poente[idx].iso.split(' ')[1][0:5]))
        output.close()
        
#####################################################################

arquivo = 'alvos.txt'				#### arquivo de alvos
horain = '2015-03-08 18:00:00'			#### hora inicial local da observacao
horafin = '2015-03-09 06:00:00'			#### hora final local da observacao
intinf = 30					#### intervalo entre as informacoes (minutos)
fuso = -3					#### fuso horario do local
latitude = '-22 32 7.8'				#### latitude do local
longitude = '314 25 2.5'			#### longitude do local
altitude = 1864					#### altitude em metros
limalt = 30					#### limite de altura para mostrar (graus)
limdist = 11					#### limite de distancia para field-of-view (arcmin)

#######################################################################

observation = Observation()
dados = observation.read(arquivo)
observation.set_sitio(longitude, latitude, altitude)
observation.set_fuse(fuso)
#observation.create_plan(horain, horafin, intinf)

