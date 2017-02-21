import numpy as np
from collections import Counter, defaultdict

adn = list('CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC')
#acomodando la cadena de adn 3'->5'
adn = adn[::-1]
#print adn
n = len(adn)
print 'ADN.length: '+ str(n)

#A1,T2,C3,G4
aminos = np.zeros(4)

for i in adn:
	if i == 'A':
		aminos[0]+=1
	elif i == 'C':
		aminos[2]+=1
	elif i == 'G':
		aminos[3]+=1
	elif i == 'T':
		aminos[1]+=1

#transition matrix
tot = len(adn)
obs = [1.0*aminos[0]/tot, 1.0*aminos[1]/tot, 1.0*aminos[2]/tot, 1.0*aminos[3]/tot]
print 'A:', obs[0], 'T', obs[1], 'C', obs[2], 'G', obs[3]
print 'A:', aminos[0], 'T', aminos[1], 'C', aminos[2], 'G', aminos[3]

aminoacidos = list('ATCG')
tm  = np.zeros((4,4))
tmCrudo = np.zeros((4,4))

desfase =Counter(zip(adn, adn[1:]))
for i in xrange(4):
	for j in xrange(4):
		t0 = aminoacidos[i]
		t1 = aminoacidos[j]
		#print t0, t1
		tm[i,j] = desfase[(t0, t1)]
		tmCrudo[i,j] = desfase[(t0, t1)]
	tm[i,] = (tm[i,])/ (sum(tm[i,]))  #lindsone smooting

print 'desfase', desfase
print 'transitionMatrix:', tm
print 'transitionMatrixCruda', tmCrudo

import random as rnd


d = {'AAAAAAAAAA': 0}
for i in xrange(100000):
	idx = np.random.choice(np.arange(4), p= obs)
	chain = aminoacidos[idx]
	for i in xrange(9):
		idx = np.random.choice(np.arange(4), p= tm[idx,])
		chain += aminoacidos[idx]
	if (chain in d):
		d[chain] += 1
		print chain
	else:
		d[chain] = 1
		print chain

import operator
#stats = {'a':1000, 'b':3000, 'c': 100, 'd':3000}
print 'Most probable chain:'
print max(d.iteritems(), key=operator.itemgetter(1))[0]