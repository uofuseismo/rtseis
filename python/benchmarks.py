#!/usr/bin/env python3
from pyrtseis.PostProcessing import Waveform as trace
from numpy import zeros
from numpy import linspace
from scipy.signal import convolve
from scipy.signal import firwin
from scipy.signal import lfilter
import time
print("begin")
trace = trace()
#print(trace.isInitialized.__doc__)
print(trace.is_initialized())
a = zeros(5000) + 1

start = time.time()
trace.set_data(a)
trace.demean()
print("Burn in time", time.time() - start, "s")

niter = 20

print("Stress it")
print("RTSeis test...");
tsum = 0
for i in range(niter):
    #start = time.time()
    trace.set_data(a+i)
    start = time.time()
    trace.demean()
    tsum = tsum + time.time() - start 
    y = trace.get_data()
    #tsum = tsum + time.time() - start
avgDemean = tsum/niter
#print("Python/C++ average time for demean   %e (s)"%(avgDemean))
tsum = 0
for i in range(niter):
   #start = time.time()
   trace.set_data(a+i)
   start = time.time()
   trace.detrend()
   tsum = tsum + time.time() - start
   y = trace.get_data()
   #tsum = tsum + time.time() - start
avgDetrend = tsum/niter

tsum = 0
for i in range(niter):
    trace.set_data(a+i)
    start = time.time()
    trace.downsample(7)
    tsum = tsum + time.time() - start
    y = trace.get_data()
avgDownsample = tsum/niter

tsum = 0
for i in range(niter):
    trace.set_data(a+i)
    start = time.time()
    trace.taper(5, "hamming")
    tsum = tsum + time.time() - start
avgTaper = tsum/niter

#print("Python/C++ average time for detrend  %e (s)"%(avgDetrend))
#trace.setData(list(1, 2, 4))
ramp = linspace(1, 5, 50)
tsum = 0
for i in range(niter):
    start = time.time()
    trace.convolve(ramp, 'full')
    tsum = tsum + time.time() - start
avgConvolve = tsum/niter
#print("Python/C++ average convolution time  %e (s)"%(avgConvolve))
#y = trace.get_data()
#print(y)
#isInitialized()

firHamming = firwin(51, 0.4, window='hamming')
tsum = 0
for i in range(niter):
    trace.set_data(a)
    start = time.time()
    trace.fir_filter(firHamming)
    tsum = tsum + time.time() - start
avgFirFilter = tsum/niter

tsum = 0
for i in range(niter):
    trace.set_data(a)
    start = time.time()
    trace.sos_lowpass_filter(0.2, order=3, prototype=trace.IIRPrototype.butterworth, zero_phase=False) #True) 
    tsum = tsum + time.time() - start
avgLPSosFilter = tsum/niter

tsum = 0
for i in range(niter):
    trace.set_data(a)
    start = time.time()
    trace.sos_bandpass_filter([0.1, 0.4], order=4, prototype=trace.IIRPrototype.butterworth, zero_phase=True)
    tsum = tsum + time.time() - start
avgBPSosFilter = tsum/niter

tsum = 0
for i in range(niter):
    trace.set_data(a+1)
    start = time.time()
    trace.envelope()
    tsum = tsum + time.time() - start
avgEnvelope = tsum/niter

print("ObsPy test...")
from obspy.core.trace import Trace
from obspy.signal.filter import envelope
traceObspy = Trace(data=a) 
traceObspy.detrend('demean')
tsum = 0
for i in range(niter):
    traceObspy.data = a + 1
    start = time.time()
    traceObspy.detrend('demean')
    tsum = tsum + time.time() - start
    y = traceObspy.data
avgDemeanObspy = tsum/niter
#print("Obspy average time for demean  %e (s)"%(avgDemeanObspy))
tsum = 0
for i in range(niter):
    traceObspy.data = a + 1 
    start = time.time()
    traceObspy.detrend('linear')
    tsum = tsum + time.time() - start
    y = traceObspy.data 
avgDetrendObspy = tsum/niter

tsum = 0
for i in range(niter):
    traceObspy.data = a + 1
    start = time.time()
    traceObspy.filter("lowpass", corners=2, freq=0.2, zerophase=False)
    tsum = tsum + time.time() - start
avgLPSosFilterObspy = tsum/niter

tsum = 0 
for i in range(niter):
    traceObspy.data = a + 1 
    start = time.time()
    traceObspy.filter("bandpass", corners=4, freqmin=0.1, freqmax=0.4, zerophase=True)
    tsum = tsum + time.time() - start
avgBPSosFilterObspy = tsum/niter 
 
#print("Obspy average time for detrend %e (s)"%(avgDetrendObspy))
tsum = 0
for i in range(niter):
    traceObspy.data = a + 1
    pct = 5./100./2.
    start = time.time()
    traceObspy.taper(pct, type='hamming')
    tsum = tsum + time.time() - start
    y = traceObspy.data
avgTaperObspy = tsum/niter

tsum = 0
for i in range(niter):
    start = time.time()
    lfilter(firHamming, 1, a)
    tsum = tsum + time.time() - start
avgFirFilterObspy = tsum/niter

tsum = 0
for i in range(niter):
    start = time.time()
    convolve(a, ramp, 'full')
    tsum = tsum + time.time() - start
avgConvolveObspy = tsum/niter

tsum = 0
for i in range(niter):
    traceObspy.data = a + 1
    start = time.time()
    traceObspy.decimate(5, no_filter=True)
    tsum = tsum + time.time() - start
avgDownsampleObspy = tsum/niter

tsum = 0
for i in range(niter):
    traceObspy.data = a + 1
    start = time.time()
    lhs = envelope(a) #traceObspy.data)
    tsum = tsum + time.time() - start
avgEnvelopeObspy = tsum/niter

#print("Average SciPy convolve time %e (s)"%(avgConvolveObspy))

print("Average demean time:     Obspy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgDemeanObspy,    avgDemean,    avgDemeanObspy/avgDemean))
print("Average detrend time:    Obspy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgDetrendObspy,   avgDetrend,   avgDetrendObspy/avgDetrend))
print("Average taper time:      Obspy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgTaperObspy,     avgTaper,     avgTaperObspy/avgTaper))
print("Average convolve time:   Scipy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgConvolveObspy,  avgConvolve,  avgConvolveObspy/avgConvolve))
print("Average FIR time:        Scipy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgFirFilterObspy, avgFirFilter, avgFirFilterObspy/avgFirFilter))
print("Average LP SOS time:     Obspy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgLPSosFilterObspy, avgLPSosFilter, avgLPSosFilterObspy/avgLPSosFilter))
print("Average BP SOS time:     Obspy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgBPSosFilterObspy, avgBPSosFilter, avgBPSosFilterObspy/avgBPSosFilter))
print("Average downsample time: Obspy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgDownsampleObspy,  avgDownsample, avgDownsampleObspy/avgDownsample))
print("Average envelope time:   Obspy %e (s), RTSeis %e (s), SpeedUp %.2fx"%(avgEnvelope, avgEnvelopeObspy, avgEnvelopeObspy/avgEnvelope))
"""
from ctypes import cdll
from ctypes import c_int
from ctypes import c_bool
from ctypes import c_float
from ctypes import c_double
from ctypes import byref
from ctypes import c_void_p
from ctypes import POINTER
from numpy import zeros
from numpy import ascontiguousarray

class Detrend:
    def __init__(self):
        lrtseis = cdll.LoadLibrary("rtseis_python.so")
        # Create the interface 
        lrtseis.rtseis_modules_detrend_initialize.argtypes = None
        lrtseis.rtseis_modules_detrend_initialize.restype = c_void_p

        lrtseis.rtseis_modules_detrend_detrend.argtypes = [c_int, POINTER(c_double), POINTER(c_double), c_void_p] 
        lrtseis.rtseis_modules_detrend_detrend.restype = c_int 
        
        lrtseis.rtseis_modules_detrend_finalize.argtypes = [c_void_p]

        # Copy the library
        self.lrtseis = lrtseis
        self.module = lrtseis.rtseis_modules_detrend_initialize()

    def __del__(self):
        self.lrtseis.rtseis_modules_detrend_finalize(self.module)
        self.module = None
        return

    def detrend(self, x):
        fname = '%s::%s'%(self.__class__.__name__, self.detrend.__name__)
        y = zeros(x.size)
        nx = len(x)
        xPtr = x.ctypes.data_as(POINTER(c_double)) 
        y = ascontiguousarray(zeros(len(x), dtype=c_double))
        yPtr = y.ctypes.data_as(POINTER(c_double))
        ierr = self.lrtseis.rtseis_modules_detrend_detrend(nx, xPtr, yPtr, self.module)
        if (ierr != 0):
            print("%s: Detrend failed")
            return None
        return y
 

from numpy import linspace
x = linspace(1, 100, 1001)
detrend = Detrend()
y = detrend.detrend(x)
"""

"""
#from numpy import zeros
#a = zeros(5)
detrend = rtseis_python.DetrendPy()
l = [1, 2, 3, 4]
detrend.detrend(l)
"""
