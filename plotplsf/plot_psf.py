import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

arq_netuno = 'ucac4_Netuno_160'
arq_triton = 'ucac4_Triton_160'

dist = 30  # distancia do centro em pixel para plot

########### ler arquivos ################

f = open(arq_netuno, 'r')
net_dados = f.readlines()
f.close
f = open(arq_triton, 'r')
tri_dados = f.readlines()
f.close

xn, yn = np.loadtxt(arq_netuno, usecols=[2,3], unpack=True)
arqsn = []
for i in net_dados:
    arqsn.append(i[350:400])
arqs = np.char.array(arqsn).strip()
xt, yt = np.loadtxt(arq_triton, usecols=[2,3], unpack=True)
arqst = []
for i in net_dados:
    arqst.append(i[350:400])
arqs = np.char.array(arqst).strip()


########### abre fits ####################

imagem = fits.open(arqs[0])
scidata = imagem[0].data
print scidata[0].shape

a, b = np.indices(scidata[0].shape)
cn = np.sqrt((b-xn)**2 + (a-yn)**2)
ct = np.sqrt((b-xt)**2 + (a-yt)**2)
dn = np.where(cn < dist)
dt = np.where(ct < dist)

########### plot ########################

f, (ax1, ax2) = plt.subplots(1, 2, sharex=True)
print len(ct[dt]), len(scidata[0][dt])

ax1.plot(cn[dn], scidata[0][dn], 'b.')
ax2.plot(ct[dt], scidata[0][dt], 'b.')
plt.savefig('plot.png', fmt='png')
