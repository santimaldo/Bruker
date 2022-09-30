# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 00:18:45 2020

@author: santi
"""

import numpy as np
import matplotlib.pyplot as plt


acq = 10
Npts = 128*2

freq = [2,3,4]
t = np.linspace(0,acq,Npts)

fid = 0
for f in freq:
   fid = fid + np.exp(-f*t) + 0.05 * (np.random.random(Npts)-0.5)
  

   
plt.figure(1)
plt.plot(t,fid)

plt.figure(2)
plt.plot(np.fft.fftshift(fid))