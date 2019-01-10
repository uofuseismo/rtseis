#!/usr/bin/env python3
from scipy import signal
from numpy import asarray
from numpy import linspace
from numpy import zeros
from math import pi
from numpy import real, imag
nw = 50
b = asarray([1611.7315706,  0.,  0.,  0.,  0.])
a = asarray([1.00000000e+00,   8.57530241e+00,   1.57767906e+02,
             7.98628595e+02,   4.76375068e+03,   7.98628595e+03,
             1.57767906e+04,   8.57530241e+03,   1.00000000e+04])

x1 =-1
x2 = 1
dx = (x2 - x1)/(nw - 1)
w = zeros(nw)
for i in range(nw):
    w[i] = 2.0*pi*pow(10, x1 + i*dx)
[win, h] = signal.freqs(b, a, w)
hr = real(h)
hi = imag(h)
for i in range(len(h)):
    print("    href1[%d] = std::complex<double> (%+.14lf,%+.14lf);"%(i, hr[i], hi[i]))

nw = 41
f = linspace(0, pi, nw)
bz = asarray([0.056340000000000, -0.000935244000000, -0.000935244000000,  0.056340000000000])
az = asarray([1.000000000000000, -2.129100000000000, 1.783386300000000, -0.543463100000000])
[fin, hz] = signal.freqz(bz, az, f)
hr = real(hz)
hi = imag(hz)
for i in range(len(fin)):
    print("    href2[%d] = std::complex<double> (%+.14lf,%+.14lf);"%(i, hr[i], hi[i]))

