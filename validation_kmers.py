from sys import argv
from collections import Counter,defaultdict
import numpy as np


window = 10	
misses = 2 
chain = list('CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC')

#Smoothing de Lindstone
def pr (x, cond, N, l=1):
	return (x+l) / (cond + l*N) 

def hammingDistance(a, b):
	n = len(a)
	count = 0
	for i in xrange(n):
		if(a[i] != b[i]):
			count += 1
	return count

def base2num(x):
    return {
        'A': 0,
        'T': 1,
        'C': 2,
        'G': 3
    }[x]	


bases = {'A':0, 'T':1, 'C':2, 'G':3}
n = len(chain)
n_bases = len(bases)

#Matriz de observaciones
obs = np.zeros(n_bases)
for i in chain:
	obs[bases[i]] += 1

nobs = np.zeros(n_bases)
for i in xrange(n_bases):
	nobs[i] = pr(obs[i], sum(obs), n_bases) 

#matriz de transiciones
tranx = np.zeros((4,4))
chains = Counter(zip(chain,chain[1:]))
for (t,t_ant), c_ws in chains.iteritems():
	tranx[bases[t], bases[t_ant]] = c_ws

#normalization
ntranx = np.zeros((4,4))
for i in xrange(4):	
	aux = sum(tranx[i,:])
	cond = n_bases
	if  (aux != 0):
		for j in xrange(4):
			ntranx[i,j] = pr(tranx[i,j], aux, cond)

print obs
print nobs
print tranx
print ntranx

########################################################
d = dict()
iteraciones = 10000
contador = 0

bases1 = list('ATCG')

while iteraciones > contador:
	prop = ['A']*window
	for j in xrange(window):
		if j == 0:
			aux = np.random.choice(4, p = nobs)
			prev = aux
			prop[0] = bases1[aux]
		else:
			aux = np.random.choice(4, p=ntranx[prev,:])
			prev = aux
			prop[j] =  bases1[aux]
		#print prop, prop[j]

	#counting number of valid k-mers
	for i in xrange(n - window+1):
		aux = chain[i:i+window]
		if hammingDistance(aux, list(prop)) <= misses:
			temp = ''.join(prop)
		 	if temp in d:
		 		d[temp] += 1
		 	else:
		 		d[temp] = 1
		 	#print prop, aux, hammingDistance(prop, aux)
	contador += 1
	#print ''


d_view = [ (v,k) for k,v in d.iteritems() ]
d_view.sort(reverse=False) # natively sort tuples by first element
for v,k in d_view:
    print "%s: %d" % (k,v)


