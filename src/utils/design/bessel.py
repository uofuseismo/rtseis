#!/usr/bin/env python3
from scipy.signal import besselap

for n in range(1,26):
   [zs,ps,k] = besselap(n)
   if (n == 1):
      print("    if (n == 1)\n    {")
   else:
      print("    }\n    else if (n == %d)\n    {"%(n))
   i = 0
   for p in ps:
      #print "        p[%d] = %.16f + %.16f*_Complex_I;"%(i,p.real,p.imag)
      print("        poles_[%d] = (%.18f, %.18f);"%(i,p.real,p.imag))
      i = i + 1
print("    }")

