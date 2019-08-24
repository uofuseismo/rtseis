#!/usr/bin/env python3
import numpy as np
import obspy

np.random.seed(4083)
vertical = np.random.randn(100)
north = np.random.randn(len(vertical))
east = np.random.randn(len(vertical))
baz = 131
radial, transverse = obspy.signal.rotate.rotate_ne_rt(north, east, baz)
aoi = 72
baz = 333
l, q, t = obspy.signal.rotate.rotate_zne_lqt(vertical, north, east, baz, aoi)

ofl = open("rotate_zne_rt_lqt.txt", "w")
for i in range(len(vertical)):
    ofl.write("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n"%( 
              vertical[i], north[i], east[i],
              radial[i], transverse[i],
              l[i], q[i], t[i]))
ofl.close()
