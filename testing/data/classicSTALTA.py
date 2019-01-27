#!/usr/bin/python
from obspy.core import read
from obspy.signal.trigger import classic_sta_lta
trace = read("https://examples.obspy.org/ev0_6.a01.gse2")[0]
df = trace.stats.sampling_rate
cft = classic_sta_lta(trace.data, int(5 * df), int(10 * df))
ofl = open('classicSTALTA_ref.txt', 'w')
for i in xrange(len(cft)):
   ofl.write("%.10f\n"%cft[i])
ofl.close()
