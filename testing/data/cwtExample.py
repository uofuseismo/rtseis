#!/usr/bin/env python3
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

if __name__ == "__main__":
    ifl = open('zwave_cwt_example.txt', 'r')
    cdat = ifl.read()
    ifl.close()
    cdat = cdat.split('\n')

    z = np.zeros(len(cdat)-1)
    for i in range(len(z)):
        z[i] = float(cdat[i])
   
    dt = 0.01
    fs = 1/dt
    fmin = 1.
    fmax = 20. # Not quite nyquist which is 50 Hz
    w0 = 6 # CWT parameter
    nf = 25
    freqs = np.linspace(fmin, fmax, nf)
    times = np.linspace(0, (len(z) - 1)*dt, len(z))

    # Scales
    widths = (w0*fs)/(2*freqs*np.pi)
    # scipy will normalize by 1/sqrt(a) = sqrt(2*f*pi*dt/w0) however.
    # obspy multiplies by convolution by dt (which it should) to get integral
    # from samples to time.  Hence, to match obspy, scipy would have to
    # multiply factor*sqrt(dt) = dt, i.e., the missing factor is factor is sqrt(dt)
    cwtmatr = signal.cwt(z, signal.morlet2, widths, w=w0)*np.sqrt(dt)

    print(cwtmatr.shape)
    ofl = open('z_cwt_amp_phase.txt', 'w')
    for j in range(len(freqs)):
        for i in range(len(z)):
            ofl.write('%f %f %e %e\n'%(times[i], freqs[j], np.abs(cwtmatr[j,i]), np.angle(cwtmatr[j,i]) ))
        ofl.write('\n')
    ofl.close()
    """
    # Rescale for plotting reasons
    cwtmatr = np.abs(cwtmatr)
    cwtmatr = (cwtmatr - np.amin(cwtmatr))/( np.amax(cwtmatr) - np.amin(cwtmatr) )
    plt.figure(figsize=(8,8))
    extent = [min(times), max(times), fmin, fmax]
    plt.imshow(cwtmatr[::-1,:], extent=extent,
               cmap='viridis', aspect='auto',
               vmin=0, vmax=1)
    pick_index = 368
    plt.plot([pick_index*dt, pick_index*dt], [-3, fmax], color='blue', linewidth=0.7)
    plt.plot(times, -2 + z*2, color='black', linewidth=0.7)
    plt.grid()
    plt.show()
    """
