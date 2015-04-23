# coding=UTF-8

import numpy as np
import itertools
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5, EarthLocation, Angle

######################################################################

ra_formatter = lambda x: "%07.4f" %x
dec_formatter = lambda x: "%06.3f" %x
alt_formatter = lambda x: "%.1f" %x
int_formatter = lambda x: "%02d" %x
int2_formatter = lambda x: "%+03d" %x

def read(datafile, name_col=None, coord_col=None, comment_col=None, time_col=None, time_fmt='jd'):
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
    if name_col != None:
        if type(name_col) not in [list, tuple, np.ndarray]:
            name_col = [name_col]
        nomes = np.loadtxt(datafile, usecols=(name_col), unpack=True, dtype ='S30', ndmin=1)
        if len(name_col) > 1:
            nome = nomes[0]
            for i in np.arange(len(name_col))[1:]:
                nome = np.core.defchararray.add(nome, ' ')
                nome = np.core.defchararray.add(nome, nomes[i])
        else:
            nome = nomes
        retornar.append(nome)
    if coord_col != None:
        if type(coord_col) not in [list, tuple, np.ndarray]:
            coord_col = [coord_col]
        coords = np.loadtxt(datafile, usecols=(coord_col), unpack=True, dtype ='S20', ndmin=1)
        coor = coords[0]
        for i in np.arange(len(coord_col))[1:]:
            coor = np.core.defchararray.add(coor, ' ')
            coor = np.core.defchararray.add(coor, coords[i])
        coord = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))
        retornar.append(coord)
    if comment_col != None:
        if type(comment_col) not in [list, tuple, np.ndarray]:
            comment_col = [comment_col]
        comments = np.loadtxt(datafile, usecols=(comment_col), unpack=True, dtype ='S30', ndmin=1)
        if len(comment_col) > 1:
            comment = comments[0]
            for i in np.arange(len(comment_col))[1:]:
                comment = np.core.defchararray.add(comment, ' ')
                comment = np.core.defchararray.add(comment, comments[i])
        else:
            comment = comments
        retornar.append(comment)
    if time_col != None:
        fmt_time = {'iso': ['S20'], 'jd': ['f8']}
        if type(time_col) not in [list, tuple, np.ndarray]:
            time_col = [time_col]
        times = np.loadtxt(datafile, usecols=(time_col), unpack=True, dtype =fmt_time[time_fmt][0], ndmin=1)
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
    """
    """
    coord = coord_pack(coord)
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
                mean_coord, dist_center = midpoint_coord(coord[a])
                if np.isscalar(dist_center.value):
                    dist_center = [dist_center]
                if all(k <= size/2 for k in dist_center):
                    same_fov.append(a)
    return same_fov

def midpoint_coord(coord, weighted=False, weight=None):
    """
    """
    i = coord_pack(coord)
    if weighted == False:
        x = (np.max(i.cartesian.x) - np.min(i.cartesian.x))/2 + np.min(i.cartesian.x)
        y = (np.max(i.cartesian.y) - np.min(i.cartesian.y))/2 + np.min(i.cartesian.y)
        z = (np.max(i.cartesian.z) - np.min(i.cartesian.z))/2 + np.min(i.cartesian.z)
    else:
        if weight == None:
            weight = [1]*len(i)
        x = np.sum(i.cartesian.x*weight)/np.sum(weight)
        y = np.sum(i.cartesian.y*weight)/np.sum(weight)
        z = np.sum(i.cartesian.z*weight)/np.sum(weight)
    delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))
    alfa = np.arctan2(y, x)
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
    ra, ts = np.meshgrid(coord.ra.deg, timeut.sidereal_time('mean').deg )
    dif_h_sid = Angle((ra-ts)*u.deg)
    dif_h_sid.wrap_at('180d', inplace=True)
    dif_h_sol = dif_h_sid * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
    dif = TimeDelta(dif_h_sol.hour*u.h, scale='tai')
    culminacao = timeut + dif
    culminacao.delta_ut1_utc = 0
    culminacao.location = site
    if rise_set == True:
        hangle_lim = np.arccos((np.cos(90.0*u.deg-limalt) - np.sin(coord.dec)*np.sin(site.latitude)) / (np.cos(coord.dec)*np.cos(site.latitude)))
        tsg_lim = Angle(ra*u.deg + hangle_lim)
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
        return culminacao, sunrise, sunset, alwaysup, neverup
    culminacao = culminacao + fuse
    return culminacao

def height_time(coord, time, time_left=False, limalt=0.0*u.deg, site=EarthLocation(0.0, 0.0, 0.0), fuse=TimeDelta(0, format='sec', scale='tai')):
    """
    """
    coord = coord_pack(coord)
    timeut = time - fuse
    if len(time.shape) == 1:
        timeut = Time([[i] for i in timeut.jd], format='jd', scale='utc')
    timeut.location = site
    timeut.delta_ut1_utc = 0
    ra, ts = np.meshgrid(coord.ra.deg, timeut.sidereal_time('mean').deg )
    hourangle = Angle((ts-ra)*u.deg)
    distzen = np.arccos(np.sin(coord.dec)*np.sin(site.latitude) + np.cos(coord.dec)*np.cos(site.latitude)*np.cos(hourangle))
    altura = 90*u.deg - distzen
    if time_left == True:
        poente = sky_time(coord, time, rise_set=True, limalt=limalt, site=site, fuse=fuse)[2]
        time_rest = poente - timeut
        return altura, time_rest
    return altura

def instant_list(time_begin, time_end=None, time_step=TimeDelta(60*60, format='sec')):
    if time_end == None:
        time_end = time_begin
    a = []
    tempo = time_begin
    while tempo <= time_end + TimeDelta(1,format='sec'):
        a.append(tempo.iso)
        tempo = tempo + time_step
    return Time(a, format='iso', scale='utc')
    
def resume_night(coord, time, name, limalt=0.0*u.deg, site=EarthLocation(0.0, 0.0, 0.0), fuse=TimeDelta(0, format='sec', scale='tai')):
    timeut = time - fuse
    coord = coord_pack(coord)
    coord_prec = precess(coord, timeut)
    culminacao, nascer, poente, alwaysup, neverup = sky_time(coord_prec, time, rise_set=True, limalt=limalt, site=site)
    ra, dec = text_coord(coord)
    night = '\nRA: ' + ra + ', DEC: ' +  dec + ', Rise: ' + np.char.array(nascer.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Culmination: ' + \
np.char.array(culminacao.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + 'TL, Set: ' + np.char.array(poente.iso).rpartition(' ')[:,:,2].rpartition(':')[:,:,0] + \
' TL, ' + np.char.array(name)
    return night
    
def text_coord(coord):
    coord = coord_pack(coord)
    ra = np.char.array([int_formatter(j) for j in coord.ra.hms.h]) + ' ' + np.char.array([int_formatter(j) for j in coord.ra.hms.m]) + ' ' + np.char.array([ra_formatter(j) for j in coord.ra.hms.s])
    sign = np.char.array(np.sign(coord.dec)).replace('-1.0', '-').replace('1.0', '+').replace('0.0', '+')
    dec = sign + np.char.array([int_formatter(j) for j in np.absolute(coord.dec.dms.d)]) + ' ' + np.char.array([int_formatter(j) for j in np.absolute(coord.dec.dms.m)]) + ' ' + np.char.array([dec_formatter(j) for j in np.absolute(coord.dec.dms.s)])
    return ra, dec

class Observation(object):
    """
    """

    def __init__(self, fuse = 0, latitude = 0.0, longitude = 0.0, height = 0.0):
        self.set_fuse(fuse)
        self.set_site(longitude, latitude, height)
        self.set_limheight(0.0)
        self.set_limdist(0.0)
        
    def set_site(self,longitude, latitude, height=0.0*u.m):
        """
        """
        self.site = EarthLocation(longitude, latitude, height)

    def set_fuse(self, fuse):
        """
        """
        self.fuse = TimeDelta(fuse*3600, format='sec', scale='tai')
        
    def set_limheight(self, height):
        self.limheight = height*u.deg
        
    def set_limdist(self, size):
        self.limdist = size*u.arcmin

    def read(self, datafile, name_col=None, coord_col=None, comment_col=None, time_col=None, time_fmt='jd'):
        a = read(datafile, name_col, coord_col, comment_col, time_col, time_fmt)
        n = 0
        if name_col != None:
            self.names = a[n]
            n = n + 1
        else:
            self.names = ['']*np.shape(a)[1]
        if coord_col != None:
            self.coords = a[n]
            n = n + 1
        if comment_col != None:
            self.comments = a[n]
            n = n + 1
        else:
            self.comments = ['']*np.shape(a)[1]
        if time_col != None:
            self.times = a[n]
            n = n + 1
        
    def close_obj(self):
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
                midcomment.append('')
            else:
                midobj.append(self.names[i][0])
                midcoord.append(self.coords[i][0])
                midcomment.append(self.comments[i][0])
        midcoord = coord_pack(midcoord)
        self.samefov = {'coords': midcoord, 'names': np.array(midobj), 'comments': np.array(midcomment), 'fov': fov}
    
#    def precess(self, coord, time):
#        """
#        """
#        return precess(coord, time)time = Time([time.iso], format='iso', scale='utc')
        
#    def sky_time(self, coord, time, rise_set=False):
#        """
#        """
#        return sky_time(self, coord, time, rise_set=False, limalt=self.limheight, site=self.site, fuse=self.fuse)
        
#    def height_time(self, coord, time, time_left=False):
#        """
#        """
#        return height_time(coord, time, time_left=False, limalt=self.limheight, site=self.site, fuse=self.fuse)
        
    def instant_list(self, time_begin, time_end=None, time_step=TimeDelta(60*60, format='sec')):
        self.instants = Time([[i] for i in instant_list(time_begin, time_end, time_step).jd], format='jd', scale='utc', location=self.site)
        self.instants.delta_ut1_utc = 0
        
    def plan(self, now=False):
        coords = self.coords
        names = self.names
        comments = self.comments
        if hasattr(self, 'samefov'):
            coords = self.samefov['coords']
            names = self.samefov['names']
            comments = self.samefov['comments']
        if not hasattr(self, 'instants') or now == True:
            instants = instant_list(Time.now() + self.fuse)
        else:
            instants = self.instants
        coord_prec = precess(coords, self.instants[0])
        culmination, lixo, lixo2, alwaysup, neverup = sky_time(coord_prec, instants[0], limalt=self.limheight, rise_set=True, site=self.site, fuse=self.fuse)
        altura, time_rest = height_time(coord_prec, instants, limalt=self.limheight, time_left=True, site=self.site, fuse=self.fuse)
        self.titles = []
        ra, dec = text_coord(coords)
        self.obs = {}
        for i in np.arange(len(instants)):
            x = np.argsort(time_rest[i].sec)
            k = np.where(altura[i,x] >= self.limheight)
            q = x[k]
            b = '\n\n---LT: {} (UT: {}), N_objects={} ----------------------------------------------------------------'.format(instants[i,0].iso.split(' ')[1][0:5], (instants[i,0] - self.fuse).iso.split(' ')[1][0:5], len(q))
            self.titles.append(b)
            if len(q) == 0:
                self.obs[instants[i,0].iso] = {'names': '', 'comments': '', 'ra': '', 'dec': '', 'height': '', 'culmination': '', 'time_left': '', 'rest': ''}
                continue
            self.obs[instants[i,0].iso] = {'names': np.char.array(names[q]), 'comments': '(' + np.char.array(comments[q]) + ')', 'ra': ra[q], 'dec': dec[q], 'height': np.char.array([alt_formatter(j) for j in altura[i,q].value]),
'culmination': np.char.array(culmination[0,q].iso).rpartition(' ')[:,2].rpartition(':')[:,0],\
'time_left': np.char.array([int_formatter(j) for j in time_rest[i,q].sec/3600.0]) + ':' + np.char.array([int_formatter(j) for j in (time_rest[i,q].sec - (time_rest[i,q].sec/3600.0).astype(int)*3600)/60]), 'rest': ['']*len(q)}
            for j in alwaysup:
                m = np.where(q == j)
                self.obs[instants[i,0].iso]['time_left'][m] = np.char.array('Always up')
            if hasattr(self, 'samefov'):
                for t in np.arange(len(q)):
                    p = q[t]
                    if len(self.samefov['fov'][p]) == 1:
                        continue
                    g = ''
                    for z in self.samefov['fov'][p]:
                        ra1, dec1 = text_coord(self.coords[z])
                        g = g + '\n\t\t' + self.names[z] + ' (' + np.char.array(self.comments[z]) + ')\n\t\t  RA: ' + ra1 + '\tDEC: ' + dec1
                    self.obs[instants[i,0].iso]['rest'][t] = g[0]
                    self.obs[instants[i,0].iso]['comments'][t] = ''

    
    def resume_night(self):
        """
        """
        self.night = resume_night(self.coords, self.instants[0], self.names, limalt=self.limheight, site=self.site, fuse=self.fuse)
            
    def create_plan(self, horain, horafin, intinf, path='.'):
        """
        """
        tempoin = Time(horain, format='iso', scale='utc', location=self.site)
        tempofin = Time(horafin, format='iso', scale='utc', location=self.site)
        intval = TimeDelta(intinf*60, format='sec')
        self.instant_list(tempoin,tempofin,intval)
        self.close_obj()
        nome = '{}/Plano_{}.dat'.format(path, tempoin.iso.split(' ')[0])
        self.plan()
        #### imprime os dados de cada objeto para cada instante ####
        output = open(nome, 'w')
        output.write('Observational Plan to the night: {}\n\n'.format(tempoin.iso.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nMinimun height: {}\nField Size: {}\n\n'.format(self.site.latitude, self.site.longitude, self.limheight, self.limdist))
        for i in np.arange(len(self.instants)):
            obs = self.obs[self.instants[i,0].iso]
            output.write(self.titles[i])
            if len(obs['names']) > 0:
                a = '\n\n' + obs['names'] + obs['comments'] + '\n\tRA: ' + obs['ra'] + '\tDEC: ' + obs['dec'] + '\n\tHeight: ' + obs['height'] + \
' deg\n\tCulmination: ' + obs['culmination'] + ' LT\n\tTime left to reach min height: ' + obs['time_left'] + obs['rest']
                for j in a:
                    output.write(j)
        self.resume_night()
        output.write('\n---Observability of the Targets----------------------------------------------------------------\n')
        for i in self.night[0]:
            output.write(i)
        output.close()
        
#####################################################################

#arquivo = 'alvos'				#### arquivo de alvos
horain = '2015-03-08 18:00:00'			#### hora inicial local da observacao
horafin = '2015-03-09 06:00:00'			#### hora final local da observacao
intinf = 60					#### intervalo entre as informacoes (minutos)
fuso = -3					#### fuso horario do local
latitude = '-22 32 7.8'				#### latitude do local
longitude = '314 25 2.5'			#### longitude do local
altitude = 1864					#### altitude em metros
#limalt = 30.0					#### limite de altura para mostrar (graus)
#limdist = 11					#### limite de distancia para field-of-view (arcmin)

########################################################################

#observation = Observation(fuso, latitude, longitude, altitude)
#objs, coords, comments = observation.read(arquivo, [0], [1,2,3,4,5,6], [7])
#observation.create_plan(coords,objs, comments, horain, horafin, intinf, limalt=limalt, size=limdist)
#print observation.plan(coord, tempoin, objs, comments, limalt=30.0*u.deg, size=11.0*u.arcmin)

a = Time.now()
obs = Observation(fuso, latitude, longitude, altitude)
obs.read('alvos', name_col=[0], coord_col=[1,2,3,4,5,6], comment_col=[7])
obs.set_limheight(30)
obs.set_limdist(11)
obs.close_obj()
##obs.instant_list(Time('2015-03-08 18:00:00', format='iso', scale='utc', location=obs.site), Time('2015-03-08 20:00:00', format='iso', scale='utc', location=obs.site))
obs.create_plan(horain, horafin, intinf)
b = Time.now()
print 'tempo de reducao: {}'.format((b-a).sec)
