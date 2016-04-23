import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import random

arq_netuno = 'ucac4_Netuno_IAG'
arq_triton = 'ucac4_Triton_IAG'

dist = 15  # distancia do centro em pixel para plot

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
arqs1 = np.char.array(arqsn).strip()
xt, yt = np.loadtxt(arq_triton, usecols=[2,3], unpack=True)
arqst = []
for i in net_dados:
    arqst.append(i[350:400])
arqs2 = np.char.array(arqst).strip()

a, b = np.indices((len(arqs1), len(arqs2)))
c = np.where(arqs1[a] == arqs2[b])


########### abre fits ####################
img = random.choice(np.arange(len(c[0])))
print img
imagem = fits.open(arqs1[c[0]][img])
scidata = imagem[0].data
print scidata.shape

a, b = np.indices(scidata.shape)
print a.shape, b.shape, xn.shape, yn.shape
cn = np.sqrt((b-xn[img])**2 + (a-yn[img])**2)
ct = np.sqrt((b-xt[img])**2 + (a-yt[img])**2)
dn = np.where(cn < dist)
dt = np.where(ct < dist)

########### plot ########################

f, (ax1, ax2) = plt.subplots(1, 2, sharex=True)
print len(ct[dt]), len(scidata[dt])

ax1.set_title("{}".format(img))

ax1.plot(cn[dn], scidata[dn], 'b.')
ax2.plot(ct[dt], scidata[dt], 'b.')
plt.savefig('plot.png', fmt='png')
