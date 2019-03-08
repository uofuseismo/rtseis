#!/bin/sh

sac << EOF
r t100.sac
taper type HAMMING width 0.2
w t100.hamming.taper.sac
r t100.sac
taper type HANNING width 0.1
w t100.hanning.taper.sac
r t100.sac
taper type COSINE width 0.3
w t100.cosine.taper.sac
r t101.sac
taper type HAMMING width 0.05
w t101.hamming.taper.sac
r t101.sac
taper type HANNING width 0.1
w t101.hanning.taper.sac
r t101.sac
taper type COSINE width 0.15
w t101.cosine.taper.sac
q
EOF 


