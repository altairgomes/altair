import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
from multiprocessing import Pool

################## Definicao de Parametros  #################################################



arquivo = 'photometry.dat'			### arquivo de dados de entrada
binsum = [1, 5, 10]					### binagem para concatenar pontos
binagemcor = [1, 3, 5, 15]				### lista de binagenspara media corrida, ps: feito apos binagem para concatenar
coltime = [2]                                   ### coluna do tempo
coltarg = [4, 7, 10, 13, 16]                    ### colunas das estrelas alvo
targcom = ['6px', '8px', '10px', '12px', '14px']			### nomes das colunas da estrela alvo - plot issue
calib = 'Yes'					### caso haja estrela de calibracao => calib = 'Yes'
colcalb = [19]					### coluna da estrela de calibracao
calbcom = ['Calib-10px']			### abertura da estrela de calibracao - plot issue
timedif = 'Yes'					### aplicar correcao de tempo
timedifv = -51.687				### correcao a ser aplicada
proc = 4					### numero maximo de processos paralelos
#  #  #  ### plot config ###  #  #  #
nomelocal = 'Maley'				### Nome do sitio de observacao
timelimuser = 'No'				### usar limites do eixo x fornecido pelo usuario => Yes ou No
timeminlim = -5					### tempo minimo do plot; segundos a partir do tempo inicial
timemaxlim = 100				### tempo maximo do plot; segundos a partir do tempo inicial
datainicio = '2013-10-25 09:40:00'		### data de referencia de contagem para o plot
fluxlimuser = 'Yes'				### usar limites de fluxo do usuario => Yes ou No
fluxminlim = 0					### fluxo minimo do plot
fluxmaxlim = 800				### fluxo maximo do plot
fluxlimcalib = 'Yes'				### usar limites de fluxo do usuario => Yes ou No
fluxmincalib = 0				### fluxo minimo do plot
fluxmaxcalib = 1200				### fluxo maximo do plot
fluxlimrelat = 'Yes'				### usar limites de fluxo do usuario => Yes ou No
fluxminrelat = 0.0				### fluxo minimo do plot
fluxmaxrelat = 2.0				### fluxo maximo do plot



#############################################################################################

colunas = np.hstack((coltime, coltarg))
coments = ['Time (jd)']
coments = np.hstack((coments, targcom))
if calib == 'Yes':
    colunas = np.hstack((colunas, colcalb))
    coments = np.hstack((coments, calbcom))
x = np.loadtxt(arquivo, usecols=(colunas), unpack=True)
t = Time(datainicio, format='iso', scale='utc')
te = t.iso.split(' ')
data = te[0]
inicio = te[1]

####### aplicando correcao de tempo ###############
if timedif == 'Yes':
    tempoantes = Time(x[0], format='jd', scale='utc')
    temponovo = tempoantes + TimeDelta(timedifv, format='sec')
    x[0] = timedifv.jd

###################### Funcoes ########################################################


###### ve a mediana de um array ####
def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]

#### concatena pontos ####
def concat(array, binsum):
    result = array[:(array.size // binsum) * binsum].reshape(-1, binsum).mean(axis=1)
    k = [array[(array.size / binsum) *binsum:].mean()]
    result = np.hstack((result, k))
    return result

#### faz a binagem corrida ####
def binamc(array, tam, binagemcor):
    return (sum(array[tam + i] for i in np.arange(binagemcor))) / binagemcor

#### salva figura ####
def salvafigura(array, coment, binsum, binagemcor, tempo, labely, ylimuser, yminlim, ymaxlim):
    axes=plt.gca()
    plt.title('{}, Data: {data}; {coment}; bin={binsum}pts; MC={binagemcor}pts'.format(nomelocal, data=data, coment=coment, binsum=binsum, binagemcor=binagemcor))
    plt.ylabel(labely)
    plt.xlabel('Tempo a partir de {inicio} (s)'.format(inicio=inicio))
    if timelimuser == 'Yes':
        axes.set_xlim([timeminlim,timemaxlim])
    if ylimuser == 'Yes':
        axes.set_ylim([yminlim,ymaxlim])
    plt.plot(tempo, array)
    plt.savefig('{}_{coment}_bin{binsum}_MC{binagemcor}.png'.format(nomelocal, coment=coment, binsum=binsum, binagemcor=binagemcor))
    plt.clf()


############################ Corpo principal #######################################################

#def plotar(binagemcor, tam):
def callfunc(binsum, binagemcor):
    global coments
    print binsum, binagemcor
####  target  ##################

    binado = np.vstack((concat(i, binsum) for i in x))

    elim = binagemcor - 1
    tam = np.arange(binado[0].size - elim)
    saida = np.vstack((binamc(i, tam, binagemcor)  for i in binado))
    tempo = (saida[0] - t.jd)*86400

### Calibracao   ###############

    if calib == 'Yes':
        h = binado[len(coltime):len(coltime) + len(coltarg)] / binado[len(coltime) + len(coltarg)]
        relat = np.vstack((binamc(i, tam, binagemcor)  for i in h))
        relatm = np.vstack((i / median(i) for i in relat))
        for k in np.arange(len(relatm)):
            texto = 'Relat_flux_{text}'.format(text=targcom[k])
            coments = np.hstack((coments,texto))
        saida = np.vstack((saida, relatm))

### Plots ##########################

    for k in np.arange(len(coltime), len(coltime) + len(coltarg)):
        salvafigura(saida[k], coments[k], binsum, binagemcor, tempo, 'Flux', fluxlimuser, fluxminlim, fluxmaxlim) 

    if calib == 'Yes':
        for k in np.arange(len(coltime) + len(coltarg), len(coltime) + len(coltarg) + len(colcalb)):
            salvafigura(saida[k], coments[k], binsum, binagemcor, tempo, 'Flux', fluxlimcalib, fluxmincalib, fluxmaxcalib) 
        for k in np.arange(len(coltime) + len(coltarg) + len(colcalb), len(saida)):
            salvafigura(saida[k], coments[k], binsum, binagemcor, tempo, 'Relative Flux', fluxlimrelat, fluxminrelat, fluxmaxrelat) 

### Save file  ###################

    strs = ["%16.8f"]
    b = ["%13.5f" for k in range(len(saida) - 1)]
    for j in b:
        strs.append(j)
    np.savetxt('photometry_{}_bin{binsum}_MC{binagemcor}.dat'.format(nomelocal, binsum=binsum, binagemcor=binagemcor), np.c_[saida.T], fmt=strs)

    f1=open('./photometry_{}_bin{binsum}_MC{binagemcor}.label'.format(nomelocal, binsum=binsum, binagemcor=binagemcor), 'w+')
    for k in np.arange(len(coments)):
        print >>f1, 'Coluna {i} = {texto}'.format(i=k+1,texto=coments[k])

#############################################################################################
######################   CHAMADA DA FUNCAO DIVIDO POR PROCESSOS   ###########################
#############################################################################################

arp, ars = [], []

for i in binsum:
    for k in binagemcor:
        arp.append(i)
        ars.append(k)

#pool = Pool(processes=proc)
map(callfunc, arp, ars)


