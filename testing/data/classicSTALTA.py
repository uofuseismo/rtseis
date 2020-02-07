#!/usr/bin/env python3
from obspy.signal.trigger import classic_sta_lta
import numpy as np

if __name__ == "__main__":
    ifl = open('gse2.txt', 'r')
    cdat = ifl.read()
    ifl.close()
    cdat = cdat.split("\n")
    trace = np.zeros(len(cdat)-1)
    for i in range(len(trace)):
        trace[i] = float(cdat[i])
    df = 200
    print(5*df, 10*df)
    cft = classic_sta_lta(trace.data, int(5 * df), int(10 * df))
    ofl = open('classicSTALTA_ref.txt', 'w')
    for i in range(len(cft)):
       ofl.write("%.10f\n"%cft[i])
    ofl.close()
