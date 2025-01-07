#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 14:02:08 2019

Transient effects: 
    - discover transient behaviour when convolving input and 
      impulse response.
    - discover that computing the output from a poduct of H and U
      results in the steady state solution.
    - visualise the transient.
    
@author: Jan Decuyper
"""
"""Alternatieve voorstellingen van de transfer functie
1. Construeer een FRF op basis van de volgende polen en nullen:
• polen: p = -0.1+5j, -0.1-5j, 3, nul: -2"""

import numpy as np
import matplotlib.pyplot as plt

# close all figures
plt.close('all')

def db(x):
    x = 20*np.log10(np.abs(x))
    return x

# 1: Construct TF

p1 = -1+2*np.pi*5*1j
p2 = -1-2*np.pi*5*1j

K = (2*np.pi*5)**2

N = 1e4
fs = 100

fvec = np.arange(0,fs,fs/N)
wvec = 2*np.pi*fvec
s = 1j*wvec

H = K/((s-p1)*(s-p2))

# 2: Compute impulse response

Hhalf = H.copy()
Hhalf[int(np.round(N/2)):] = 0

h = 2*np.real(np.fft.ifft(Hhalf))
tvec = np.arange(0,N/fs,1/fs)

plt.figure(1)
plt.plot(tvec,h*fs)
plt.xlabel('Time (s)')
plt.ylabel('h(t)')
plt.title('Impulse response')

# 3: Apply input sine

f = 0.5

tvec = np.arange(0,10/f,1/fs)
N = len(tvec)

u = np.sin(2*np.pi*f*tvec)
y = np.convolve(h,u)
y = y[0:int(N)]

NPP = int((1/f)*fs)

plt.figure(2)
plt.plot(tvec,u,'b')
plt.plot(tvec,y,'r')
plt.plot(tvec[N-NPP:],y[N-NPP:],'g*')
plt.xlabel('Time (s)')
plt.ylabel('y(t)')
plt.legend(('input','output'))

# 4: Show transient

LastP = y[N-NPP:] # (N-1)-NPP+1 om niet begin en einde interval mee te nemen

Steady_state = np.tile(LastP,10)
Transient = y-Steady_state

plt.figure(3)
plt.plot(tvec,y,'r')
plt.plot(tvec,Steady_state,'g:')
plt.plot(tvec,Transient,'k')
plt.xlabel('Time (s)')
plt.ylabel('y(t)')
plt.legend(('Output','Steady state','Transient'))

#%%
# z=np.array([1,2,3])
# sum1 =0
# for i in z:
#     sum1 = sum1+i
    
# spring_fruits = ['Apricot', 'Avocado', 'Kiwi', 'Grapefruit', 'Cherry', 'Strawberry']

# summer_fruits = spring_fruits.copy()

# print(summer_fruits)

