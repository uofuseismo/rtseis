#!/usr/bin/env python3
from numpy import ones
from obspy.core.trace import Trace
t100 = Trace(ones(100))
t100.write('t100.sac', format='sac')
t101 = Trace(ones(101))
t101.write('t101.sac', format='sac')
