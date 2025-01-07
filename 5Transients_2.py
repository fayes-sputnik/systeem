#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:00:22 2019
Compare the computed output in the time-domain to the output computed
in the frequency domain.

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

# close all figures
plt.close('all')

def db(x):
    x = 20*np.log10(np.abs(x))
    return x

# 1: Compute FRF, tvec, h(tvec)

# system  
p1 = -4+100*1j
p2 = -4-100*1j
z = -100
K = 100

# fvec
fs = 1000
N = 2560
fvec = np.arange(0,fs,fs/N)
wvec = 2*np.pi*fvec

# H(s)
s = 1j*wvec
H = K*(s-z)/((s-p1)*(s-p2))

# h(tvec)
Hhalf = H.copy()
Hhalf[int(np.round(N/2)):] = 0
h = 2*np.real(np.fft.ifft(Hhalf))

tvec = np.arange(0,np.around(N*(1/fs),5),1/fs)

# 2: apply input signal
T = N*(1/fs)/8
f = 1/T
u = np.cos(2*np.pi*f*tvec)

y = np.convolve(h,u)
y = y[0:N]

NPP = int(T*fs)
LastP = y[N-NPP:]
Steady_state = np.concatenate([LastP,LastP,LastP,LastP,LastP,LastP,LastP,LastP]) # last part is 8 times due to 8 periods
Transient = y-Steady_state

plt.figure(10)
plt.plot(tvec,h,'c')

plt.figure(1)
plt.plot(tvec,y,'r')
plt.plot(tvec,Steady_state,'g:')
plt.xlabel('Time (s)')
plt.ylabel('y(t)')
plt.legend(('Output','Steady state'))

plt.figure(2)
plt.plot(tvec,Transient,'k')
plt.xlabel('Time (s)')
plt.ylabel('Transient(t)')

Y = np.fft.fft(y) # unscaled spectrum
STEADY_STATE = np.fft.fft(Steady_state) # unscaled spectrum

plt.figure(3)
plt.plot(fvec,db(Y),'r+')
plt.plot(fvec,db(STEADY_STATE),'g.')
plt.xlim(0,fs/2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('|Y| (dB)')
plt.legend(('Spectrum of the output','Spectrum of the steady state'))

# 3: Frequency domain output
U = np.fft.fft(u) # unscaled spectrum
Yfreq = H*U # unscaled output
Yfreq_half = Yfreq.copy()
Yfreq_half[int(np.round(N/2)):] = 0

yfreq = 2*np.real(np.fft.ifft(Yfreq_half))

plt.figure(4)
plt.plot(fvec,db(Y),'r+')
plt.plot(fvec,db(Yfreq),'b.')
plt.xlabel('Frequency (Hz)')
plt.ylabel('|Y| (dB)')
plt.legend(('Output from convolution','Output from Y(s) = H(s)*U(s)'))
plt.xlim(0,fs/2)


plt.figure(5)
plt.plot(tvec,y,'r')
plt.plot(tvec,yfreq,'b')
plt.xlabel('Time (s)')
plt.ylabel('y(t)')
plt.legend(('y(t) from convolution','y(t) from ifft of Y(s) = H(s)*U(s)'))
plt.show()