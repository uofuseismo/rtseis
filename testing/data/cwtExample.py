#!/usr/bin/env python3
import numpy as np
from scipy import signal
#from sig_work import cwt
import matplotlib.pyplot as plt

def cwt_custom_scipy(data, wavelet, widths, dtype=None, **kwargs):
    """
    Custom implementation of scipy's CWT.  This agrees very closely with ObsPy's
    (after changing the way ObsPy chooses scales).
    """
    if wavelet == signal.ricker:
        window_size = kwargs.pop('window_size', None)
    # Determine output type
    if dtype is None:
        if np.asarray(wavelet(1, widths[0], **kwargs)).dtype.char in 'FDG':
            dtype = np.complex128
        else:
            dtype = np.float64

    output = np.zeros((len(widths), len(data)), dtype=dtype)
    for ind, width in enumerate(widths):
        N = np.min([10 * width, len(data)])
        # the conditional block below and the window_size
        # kwarg pop above may be removed eventually; these
        # are shims for 32-bit arch + NumPy <= 1.14.5 to
        # address gh-11095
        if wavelet == signal.ricker and window_size is None:
            ceil = np.ceil(N)
            if ceil != N:
                N = int(N)
        N = len(data) #Overwriting this line
        wavelet_data = np.conj(wavelet(N, width, **kwargs)[::-1])
        output[ind] = np.convolve(data, wavelet_data, mode='same')
    return output

if __name__ == "__main__":
    ifl = open('zwave_cwt_example.txt', 'r')
    cdat = ifl.read()
    ifl.close()
    cdat = cdat.split('\n')

    z = np.zeros(len(cdat)-1)
    for i in range(len(z)):
        z[i] = float(cdat[i])
    #z[:] = 0
    #z[0] = 1
   
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
    #cwtmatr = signal.cwt(z, signal.morlet2, widths, w=w0)*np.sqrt(dt)
    cwtmatr = cwt_custom_scipy(z, signal.morlet2, widths, w=w0)*np.sqrt(dt)

    print(cwtmatr.shape)
    ofl = open('z_cwt_amp_phase.txt', 'w')
    for j in range(len(freqs)):
        for i in range(len(z)):
            ofl.write('%f %f %e %e\n'%(times[i], freqs[j], np.abs(cwtmatr[j,i]), np.angle(cwtmatr[j,i])) )#, np.real(cwtmatr[j,i]), np.imag(cwtmatr[j,i]) ))
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
