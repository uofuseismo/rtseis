#!/usr/bin/env python3
import numpy as np
from obspy.signal.interpolation import weighted_average_slopes
import matplotlib.pyplot as plt
# read data
old_start = 0
old_dt = 1/200.
new_start = 0
new_dt = 1./250.0
ifl = open('gse2.txt', 'r')
cdat = ifl.read()
cdat = cdat.split('\n')
ifl.close()
# unpack data
n = len(cdat) - 1
x = np.zeros(n)
for i in range(n):
    x[i] = float(cdat[i])
# interpolate
tmax = (len(x) - 1)*old_dt
new_npts = int(tmax/new_dt + 0.5)
while (new_npts - 1)*new_dt > tmax:
    new_npts = new_npts - 1
xint = weighted_average_slopes(x, old_start, old_dt, new_start, new_dt, new_npts)
ofl = open('wigint.txt', 'w')
for i in range(len(xint)):
    ofl.write('%f, %.8f\n'%(i*new_dt, xint[i]))
ofl.close()
#print(xint)
lplot = False
if (lplot):
    t_old = np.linspace(0, (n - 1)*old_dt, n) 
    t_new = np.linspace(0, (new_npts - 1)*new_dt, new_npts)
    plt.plot(t_old, x, color='black')
    plt.plot(t_new, xint, linestyle='--', color='red', linewidth=0.75)
    plt.xlim(33,35)
plt.show()
