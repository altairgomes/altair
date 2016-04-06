import numpy as np

a = 'ucac4_Triton_160'

f = open(a, 'r')
dados = f.readlines()
f.close()

filtro, exp = [], []

limexp = 10

filtros_zei = {'': 'desconhecido', 'Clear': 'clear', 'I': 'I', 'R': 'R', 'filtro desconhecido': 'desconhecido'}
filtros_iag = {'': 'desconhecido', "'R       '": 'R', 'B': 'B', 'CLEAR': 'Clear', 'Clear': 'Clear', 'Dark': 'desconhecido', 'I': 'I', 'Metano': 'Metano', 'R': 'R',
               'U': 'U', 'V': 'V', 'clear': 'Clear', 'filtro desconhecido': 'desconhecido'}
filtros_160 = {'': 'desconhecido', '-': 'desconhecido', '14/07/01': 'desconhecido', '15/07/01': 'desconhecido', '5magn': 'desconhecido', 'B': 'B', 'C': 'Clear', 'CLEAR': 'Clear', 'Cl': 'Clear',
               'Clear': 'Clear', 'CuSO4': 'desconhecido', 'I': 'I', 'Metano': 'Metano', 'R': 'R', 'SEM': 'desconhecido', 'V': 'V', 'cl': 'Clear', 'clear': 'Clear','filtro desconhecido': 'desconhecido'}

filtros = filtros_160

for i in dados:
    filtro.append(i[327:349])
    exp.append(i[321:327])
    
filtro = np.char.array(filtro).strip()
exp = np.array(exp, dtype=np.int)

#print np.unique(filtro)

for i in filtros.keys():
    n = np.where((filtro == i) & (exp > limexp))
    f = open(a + '_' + filtros[i] + '_' + 'b10', 'a')
    for k in n[0]:
        f.write(dados[k])
    f.close()
    n = np.where((filtro == i) & (exp <= limexp))
    f = open(a + '_' + filtros[i] + '_' + 's10', 'a')
    for k in n[0]:
        f.write(dados[k])
    f.close()

