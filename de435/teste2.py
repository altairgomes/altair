import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
import itertools
import matplotlib.pyplot as plt

a = ['Ananke', 'Carme', 'Elara', 'Himalia', 'Leda', 'Lysithea', 'Pasiphae', 'Sinope']
b = [874, 511, 809]
c = ['PE', 'BC', 'OH', 'E', 'Z']

plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]

params = {'text.usetex' : True,
          'font.size' : 22,
          'font.family' : 'txfonts',
          'text.latex.unicode': True,
          }
          
plt.rcParams.update(params)

sizel = 18
sizek = 14

def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        >>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    data.sort()
    groups = [[0]]
    for i in np.arange(len(data[1:]))+1:
        if abs(data[1:][i-1] - data[groups[-1][-1]]) <= maxgap:
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups
    
    
for i in a:
    print i
    for k in b:
        data = np.loadtxt('{}_{}_off.dat'.format(i,k), usecols=[0,1,2,3,4], dtype={'names': ('ra','dec','jd','av','tel'), 'formats': ('f8', 'f8', 'f16', 'f8','S2')})
        for l in c:
            tel = np.where(data['tel'] == l)
            if len(tel[0]) == 0:
                continue
            g = cluster(data['jd'][tel[0]], TimeDelta(10*u.h).jd)
            ra = np.array([data['ra'][tel[0]][x].mean() for x in g])
            era = np.array([data['ra'][tel[0]][x].std() for x in g])
            av = np.array([data['av'][tel[0]][x].mean() for x in g])
            if l == 'PE':
                plt.errorbar(av, ra*1000, yerr=era*1000, fmt='s', label='PE', color='red')
            if l == 'BC':
                plt.errorbar(av, ra*1000, yerr=era*1000, fmt='o', label='B&C', color='blue')
            if l == 'OH':
                plt.errorbar(av, ra*1000, yerr=era*1000, fmt='^', label='OHP', color='black')
            if l == 'E':
                plt.errorbar(av, ra*1000, yerr=era*1000, fmt='*', label='ESO', color='green')
            if l == 'Z':
                plt.errorbar(av, ra*1000, yerr=era*1000, fmt='v', label='Zeiss', color='magenta')
    plt.title('{}: Right Ascension, DE435'.format(i), fontsize=sizel)
    plt.xlim(0,360)
    plt.ylim(-600,600)
    plt.xlabel('True Anomaly', fontsize=sizel)
    plt.xticks(np.arange(0,361,45))
    plt.ylabel('Offset (mas)', fontsize=sizel)
    plt.legend(numpoints=1, prop={'size':sizek})
    plt.axhline(0, color='black')
    fig = plt.gcf()
    fig.set_size_inches(7.2,4.05)
    plt.tick_params(axis='both', which='major', labelsize=sizel)
    fig.savefig('{}_RA.png'.format(i),dpi=300, format='png', bbox_inches='tight')
    plt.clf()
    for k in b:
        data = np.loadtxt('{}_{}_off.dat'.format(i,k), usecols=[0,1,2,3,4], dtype={'names': ('ra','dec','jd','av','tel'), 'formats': ('f8', 'f8', 'f16', 'f8','S2')})
        for l in c:
            tel = np.where(data['tel'] == l)
            if len(tel[0]) == 0:
                continue
            g = cluster(data['jd'][tel[0]], TimeDelta(10*u.h).jd)
            dec = np.array([data['dec'][tel[0]][x].mean() for x in g])
            edec = np.array([data['dec'][tel[0]][x].std() for x in g])
            av = np.array([data['av'][tel[0]][x].mean() for x in g])
            if l == 'PE':
                plt.errorbar(av, dec*1000, yerr=edec*1000, fmt='s', label='PE', color='red')
            if l == 'BC':
                plt.errorbar(av, dec*1000, yerr=edec*1000, fmt='o', label='B&C', color='blue')
            if l == 'OH':
                plt.errorbar(av, dec*1000, yerr=edec*1000, fmt='^', label='OHP', color='black')
            if l == 'E':
                plt.errorbar(av, dec*1000, yerr=edec*1000, fmt='*', label='ESO', color='green')
            if l == 'Z':
                plt.errorbar(av, dec*1000, yerr=edec*1000, fmt='v', label='Zeiss', color='magenta')
    plt.title('{}: Declination, DE435'.format(i), fontsize=sizel)
    plt.xlim(0,360)
    plt.ylim(-600,600)
    plt.xlabel('True Anomaly', fontsize=sizel)
    plt.xticks(np.arange(0,361,45))
    plt.ylabel('Offset (mas)', fontsize=sizel)
    plt.legend(numpoints=1, prop={'size':sizek})
    plt.axhline(0, color='black')
    fig =plt.gcf()
    fig.set_size_inches(7.2,4.05)
    plt.tick_params(axis='both', which='major', labelsize=sizel)
    fig.savefig('{}_DEC.png'.format(i),dpi=300, format='png', bbox_inches='tight')
    plt.clf()
