#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import CubicHermiteSpline

def weighted_average_slopes(data, old_start, old_dt, new_start, new_dt,
                            new_npts):#, *args, **kwargs):
    # In almost all cases the unit will be in time.
    new_end = new_start + (new_npts - 1)*new_dt
    new_time_array = np.linspace(new_start, new_end, new_npts)

    m = np.diff(data) / old_dt
    w = np.abs(m)
    w = 1.0 / np.clip(w, np.spacing(1), w.max())

    slope = np.empty(len(data), dtype=np.float64)
    slope[0] = m[0]
    slope[1:-1] = (w[:-1] * m[:-1] + w[1:] * m[1:]) / (w[:-1] + w[1:])
    slope[-1] = m[-1]

    # If m_i and m_{i+1} have opposite signs then set the slope to zero.
    # This forces the curve to have extrema at the sample points and not
    # in-between.
    sign_change = np.diff(np.sign(m)).astype(np.bool)
    slope[1:-1][sign_change] = 0.0

    derivatives = np.empty((len(data), 2), dtype=np.float64)
    print(data[0:30])
    print(slope[0:30])
    derivatives[:, 0] = data
    derivatives[:, 1] = slope
    x = np.linspace(old_start, (len(x) - 1)*old_dt, len(x))
    interp = CubicHermiteSpline(x, data, derivatives)

def load_data(fname = 'gse2.txt'):
    ifl = open(fname, 'r')
    cdat = ifl.read()
    cdat = cdat.split('\n')
    ifl.close()
    # unpack data
    n = len(cdat) - 1
    x = np.zeros(n)
    for i in range(n):
        x[i] = float(cdat[i]) 
    return x

"""
old_start = 0
old_dt = 1/200.
new_start = 0
new_dt = 1./250.0
x = load_data('gse2.txt')
nq_new = 14999
weighted_average_slopes(x, old_start, old_dt, new_start, new_dt, nq_new)
"""
                         
from obspy.signal.interpolation import weighted_average_slopes
import matplotlib.pyplot as plt
# read data
old_start = 0
old_dt = 1/200.
new_start = 0
new_dt = 1./250.0
x = load_data('gse2.txt')
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
