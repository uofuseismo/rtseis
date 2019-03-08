#!/usr/bin/env python3
from obspy.core import read
from numpy import sqrt
t1_hamming = read('t100.hamming.taper.sac')
t1_hanning = read('t100.hanning.taper.sac')
t1_cosine  = read('t100.cosine.taper.sac')

t2_hamming = read('t101.hamming.taper.sac')
t2_hanning = read('t101.hanning.taper.sac')
t2_cosine  = read('t101.cosine.taper.sac')

cline1 = ''
for i in range(len(t1_hamming[0].data)):
    cline1 = cline1 + '%.8lf, %.8lf, %.8lf\n'%(t1_hamming[0].data[i], t1_hanning[0].data[i], pow(t1_cosine[0].data[i],1))
print(cline1)
of = open('../taper100.all.txt', 'w')
of.write(cline1)
of.close()

cline2 = ''
for i in range(len(t2_hamming[0].data)):
    cline2 = cline2 + '%.8lf, %.8lf, %.8lf\n'%(t2_hamming[0].data[i], t2_hanning[0].data[i], pow(t2_cosine[0].data[i],1))
print(cline2)
of = open('../taper101.all.txt', 'w')
of.write(cline2)
of.close()


"""
of = open('../taper100.hanning.txt', 'w')
of.write(cline1)
of.close()

cline2 = ''
for d in t2[0].data:
    cline2 = cline2 + '%lf\n'%d
print(cline2)
of = open('../taper101.hamming.txt', 'w')
of.write(cline2)
of.close()
"""
