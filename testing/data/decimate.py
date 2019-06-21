#!/usr/bin/env python3
from scipy import signal
import numpy as np

def load_data(filename = '../data/gse2.txt'):
    ifl = open(filename, 'r')
    cdat = ifl.read()
    ifl.close()
    cdat = cdat.split("\n")
    data = np.zeros(len(cdat)-1)
    for i in range(len(data)):
        data[i] = float(cdat[i])
    return data

def fir_decimate_and_write(data, nq, filter_len, zero_phase=True):
    dec = signal.decimate(data, nq, n=filter_len-1,
                          ftype='fir', zero_phase=zero_phase) 
    ofl_name = 'decimate_%d.txt'%nq
    ofl = open(ofl_name, 'w')
    for i in range(len(dec)):
        ofl.write('%.12f\n'%dec[i])
    ofl.close() 

if __name__ == "__main__":
    data = load_data()
    fir_decimate_and_write(data, 2,  97)
    fir_decimate_and_write(data, 3,  97) 
    fir_decimate_and_write(data, 4,  97) 
    fir_decimate_and_write(data, 5,  101) 
    fir_decimate_and_write(data, 6,  97) 
    fir_decimate_and_write(data, 7,  99) 
    fir_decimate_and_write(data, 8,  97) 
    fir_decimate_and_write(data, 9,  109) 
    fir_decimate_and_write(data,10,  101) 







 
"""
sc.signal.decimaet( ) 
function decimateTest()
    % Load the reference data
    data = load("../data/gse2.txt");
    d2 = decimate(data, 2, 97-1, 'FIR');
    d3 = decimate(data, 3, 109+1, 'FIR');
    d4 = decimate(data, 4, 97+1,  'FIR');
    d5 = decimate(data, 5, 99+1,  'FIR');
    d6 = decimate(data, 6, 97,  'FIR');
    d7 = decimate(data, 7, 101, 'FIR');
    d8 = decimate(data, 8, 97,  'FIR');
    d9 = decimate(data, 9, 97,  'FIR');
    d10= decimate(data, 10, 97, 'FIR');
"""
