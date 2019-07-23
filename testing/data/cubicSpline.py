#!/usr/bin/env python3
from scipy.interpolate import CubicSpline
import numpy as np

x = np.arange(21)*np.pi/5.0
y = np.sin(x)
y[0] = 0.01
y[-1] = y[0]
xq = np.linspace(min(x), max(x), 101)
cs_not_a_knot = CubicSpline(x, y, bc_type='not-a-knot') 
cs_clamped    = CubicSpline(x, y, bc_type='clamped')
cs_natural    = CubicSpline(x, y, bc_type='natural')
cs_periodic   = CubicSpline(x, y, bc_type='periodic')
yq_knot    = cs_not_a_knot(xq)
yq_clamped = cs_clamped(xq)
yq_natural = cs_natural(xq)
yq_period  = cs_periodic(xq)
ofl = open('cubicSplineInterpolation.txt', 'w')
for i in range(len(xq)):
    ofl.write("%e, %e, %e, %e, %e\n"%(xq[i], yq_knot[i], yq_clamped[i],
                                      yq_natural[i], yq_period[i]))
ofl.close()
print(cs_not_a_knot.integrate(0, 1))
print(cs_natural.integrate(max(x), 0))
print(cs_clamped.integrate(4.08, 6))
print(cs_periodic.integrate(3.99, max(x)))
print(cs_periodic.integrate(3.007729,2.541008))
