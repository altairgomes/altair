import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob
from astropy.time import Time
import scipy.optimize as optimization

#########################################################

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
    if len(coord_col) > 0:
        n = len(coord_col)
        coords = np.arange(0, n/2, dtype=np.int8) + len(name_col)
        ra = dados[coords[0]]
        dec = dados[coords[0]+n/2]
        for i in coords[1:]:
            ra = np.core.defchararray.add(ra, ' ')
            ra = np.core.defchararray.add(ra, dados[i])
            dec = np.core.defchararray.add(dec, ' ')
            dec = np.core.defchararray.add(dec, dados[i+n/2])
        retornar['coords'] = np.array([ra, dec])#coor
    if len(time_col) > 0:
        n = len(name_col) + len(coord_col) + len(comment_col)
        times = np.arange(0, len(time_col), dtype=np.int8) + n
        if time_fmt == 'iso':  ####### iso nao deve funcionar por enquanto
            tim = dados[times[0]]
            a = len(time_col)
            len_iso = {2: [' '], 6: ['-', '-', ' ', ':',':']}
            for i in times[1:]:
                print a, i
                tim = np.core.defchararray.add(tim, len_iso[a][i-n-1]) 
                tim = np.core.defchararray.add(tim, dados[i])
        elif time_fmt == 'jd':
            tim = dados[times[0]].astype(np.float)
        retornar['times'] = tim
    return retornar

def readxy(f):
    for i in ['ucac4.red.xy', '2mass.red.xy']:
        arq = f[:-4]+i
        x, y, ra, dec, jd, s, h = np.loadtxt(arq, usecols=[0,1,33,34,41, 29, 3], unpack=True)
        if (s[0] > 2 and len(x) > 2):
            coord = SkyCoord(ra,dec, unit=(u.hourangle,u.deg))
            time = Time(jd[0], format='jd', scale='utc')
            d = np.where(np.absolute(ephnet['time'].jd - time.jd)*u.day < 0.1*u.s)
            e = np.where(coord.separation(ephnet['coord'][d]) < 1.5*u.arcsec)
            if len(e) > 0:
                ra = np.delete(ra, e)
                dec = np.delete(dec, e)
                x = np.delete(x, e)
                y = np.delete(y, e)
            d = np.where(np.absolute(ephtri['time'].jd - time.jd)*u.day < 0.1*u.s)
            e = np.where(coord.separation(ephtri['coord'][d]) < 1.5*u.arcsec)
            if len(e) > 0:
                ra = np.delete(ra, e)
                dec = np.delete(dec, e)
                x = np.delete(x, e)
                y = np.delete(y, e)
            coord = SkyCoord(ra, dec, unit=(u.hourangle,u.deg))
            return {'x': x, 'y': y, 'coord': coord, 'jd': time}
    return None
    
def proj(coord, coord0):
    a = np.sin(coord.dec)*np.sin(coord0.dec)+np.cos(coord.dec)*np.cos(coord0.dec)*np.cos(coord.ra-coord0.ra)
    X = np.cos(coord.dec)*np.sin(coord.ra - coord0.ra)/a
    Y = (np.sin(coord.dec)*np.cos(coord0.dec)-np.cos(coord.dec)*np.sin(coord0.dec)*np.cos(coord.ra-coord0.ra))/a
    return X,Y
    
def func(x, a, b, c):
    return a + b*x[0] + c*x[1]
    
def cria_mes(i, ajx, ajy, seeing):
    ### estrelas
    X, Y = proj(meancoord, dados[i]['coord'][0])
    x = (ajy[2]*X-ajx[2]*Y-ajy[2]*ajx[0]+ajx[2]*ajy[0])/(ajy[2]*ajx[1]-ajx[2]*ajy[1])
    y = (ajy[1]*X-ajx[1]*Y-ajy[1]*ajx[0]+ajx[1]*ajy[0])/(ajx[2]*ajy[1]-ajy[2]*ajx[1])
    ### netuno
    d = np.where(np.absolute(ephnet['time'].jd - dados[i]['jd'].jd)*u.day < 0.1*u.s)
    Xn, Yn = proj(ephnet['coord'][d], dados[i]['coord'][0])
    xn = (ajy[2]*Xn-ajx[2]*Yn-ajy[2]*ajx[0]+ajx[2]*ajy[0])/(ajy[2]*ajx[1]-ajx[2]*ajy[1])
    yn = (ajy[1]*Xn-ajx[1]*Yn-ajy[1]*ajx[0]+ajx[1]*ajy[0])/(ajx[2]*ajy[1]-ajy[2]*ajx[1])
    ### tritao
    d = np.where(np.absolute(ephtri['time'].jd - dados[i]['jd'].jd)*u.day < 0.1*u.s)
    Xt, Yt = proj(ephtri['coord'][d], dados[i]['coord'][0])
    xt = (ajy[2]*Xt-ajx[2]*Yt-ajy[2]*ajx[0]+ajx[2]*ajy[0])/(ajy[2]*ajx[1]-ajx[2]*ajy[1])
    yt = (ajy[1]*Xt-ajx[1]*Yt-ajy[1]*ajx[0]+ajx[1]*ajy[0])/(ajx[2]*ajy[1]-ajy[2]*ajx[1])
    ### cria arquivo mes
    f = open(i[:-4]+'mes', 'w')
    for k in np.arange(len(x)):
        if (x[k] > seeing + 1 and y[k] > seeing + 1 and x[k] < 2048 - seeing -1 and y[k] < 2048 - seeing - 1):
            f.write('circle(  {},   {},    {})\n'.format(x[k],y[k],seeing))
#    f.write('circle(  {},   {},    {})\n'.format(xn[0],yn[0],seeing))
#    f.write('circle(  {},   {},    {})'.format(xt[0],yt[0],seeing))
    f.close()    
    return x, y
    
    
#########################################################
t = Time.now()
eph = read2('../target_Netuno_160', coord_col=[0,1,2,3,4,5], time_col=[6], time_fmt='jd')
ephnet = {'coord': SkyCoord(eph['coords'][0], eph['coords'][1], unit=(u.hourangle, u.deg)), 'time': Time(eph['times'], format='jd', scale='utc')}
print 'read ephnet: ', (Time.now()-t).sec

t = Time.now()
eph2 = read2('../target_Triton_160', coord_col=[0,1,2,3,4,5], time_col=[6], time_fmt='jd')
ephtri = {'coord': SkyCoord(eph2['coords'][0], eph2['coords'][1], unit=(u.hourangle, u.deg)), 'time': Time(eph2['times'], format='jd', scale='utc')}
print 'read ephtri: ', (Time.now()-t).sec
#########################################################

lista = glob.glob('Netuno_0*.fits')

dados = {}
t = Time.now()
for i in lista:
    k = readxy(i)
    if k:
        dados.update({i: k})
print 'dados: ', (Time.now()-t).sec

meancoordra = dados[lista[0]]['coord'].ra.value
meancoorddec = dados[lista[0]]['coord'].dec.value

t = Time.now()
for i in dados.keys():
    c = SkyCoord(meancoordra, meancoorddec, unit=(u.deg,u.deg))
    d = dados[i]['coord']
    n,m = np.indices((len(c), len(d)))
    k = np.where(c[n].separation(d[m]) < 1.5*u.arcsec)
    l = [x for x in np.arange(len(d)) if x not in k[1]]
    if len(l) == 0:
        continue
    meancoordra = np.concatenate((meancoordra, d[l].ra.value))
    meancoorddec = np.concatenate((meancoorddec, d[l].dec.value))
    
print 'calcula stars: ', (Time.now()-t).sec
#print meancoordra, meancoorddec, len(meancoordra)
meancoord = SkyCoord(meancoordra, meancoorddec, unit=(u.deg,u.deg))
    
## calcular ajuste e criar arquivos mes

t = Time.now()
for i in dados.keys():
    print i
    d = dados[i]
    X,Y = proj(dados[i]['coord'], dados[i]['coord'][0])
######## otimizacao
    ajx = optimization.curve_fit(func, np.vstack((dados[i]['x'],dados[i]['y'])), X, np.zeros(3))
    ajy = optimization.curve_fit(func, np.vstack((dados[i]['x'],dados[i]['y'])), Y, np.zeros(3))
    
    pscale = d['coord'][1].separation(d['coord'][0])/(np.sqrt((d['x'][1]-d['x'][0])*(d['x'][1]-d['x'][0]) + (d['y'][1]-d['y'][0])*(d['y'][1]-d['y'][0])))
    seeing = 2.0/(pscale.to(u.arcsec).value)
    cria_mes(i, ajx[0], ajy[0], seeing)
    
print 'cria mes: ', (Time.now()-t).sec
