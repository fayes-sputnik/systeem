#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 06:32:38 2019

Correct usage of the DFT. Discovering leakage. When and how to avoid it.

@author: Jan Decuyper
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
# close all figures
plt.close('all')

"""Maak een sinus met f = 2 Hz, en bemonster deze met fs = 10 Hz. Neem N = 10, 11, · · · , 15.
Plot telkens het tijdssignaal en het frequentiespectrum. Herschaal het frequentiespectrum
met N om de invloed van het aantal punten op de schaling van het spectrum te ontwijken.
Hoe be¨ınvloed de meettijd (N) het spectrum van het signaal?"""

def db(x):
    x = 20*np.log10(np.abs(x))
    return x

# 1: sinus 2 Hz range of N

f = 2 # Hz
fs = 10.0 # Hz
Ts = 0.1 # s

# Continue sinus
fs_cont = 100 # Hz
Ts_cont = 1/fs_cont # s
tvec_cont = np.arange(-1,3,Ts_cont)
y_cont = np.sin(2*np.pi*f*tvec_cont)

N = np.array([10,11,12,13,14,15])

for i in range(len(N)):
    Ni = N[i]
    Tend = Ni*Ts
    
    # tijdsdomein
    tvec = np.arange(0,np.around(Tend,3),Ts) # round to a number of decimals to avoid machine precission errors
    y = np.sin(2*np.pi*f*tvec)
    
    plt.figure(1)
    plt.subplot(len(N),2,(2*i)+1)
    plt.plot(tvec,y,'r*-')
    plt.plot(tvec_cont,y_cont,':k')
    plt.ylabel('y(t)')
    plt.xlabel('Time (s)')
    

    # frequentiedomein
    fvec = np.arange(0,fs,fs/Ni)
    Y = np.fft.fft(y)/Ni # altijd schalen met N
    
    plt.figure(1)
    plt.subplot(len(N),2,(2*i)+2)
    plt.stem(fvec,np.abs(Y))
    plt.ylabel('|Y|')
    plt.xlabel('Frequency (Hz)')
    # 
# sys.exit()
# 2: sinus 80 Hz
"""Maak een signaal u(t) = Asin(wt + phi), n = 0, 1, · · · ,N ;bepaal DFT met schaling 1/N"""

f = 80
A = 1
fs = 1000
N = 16
phi = np.pi/2

Tend = np.around(N*1/fs,5)
tvec = np.arange(0,Tend,1/fs)

u = np.sin(2*np.pi*f*tvec+phi)    
U = np.fft.fft(u)/N
#plot u(t)
plt.figure(2)
plt.subplot(221)
plt.plot(tvec,u,'k*:')
plt.xlabel('Time (s)')
plt.ylabel('u(t)')

#plot |U(f)|
fvec = np.arange(0,fs,fs/N)

plt.figure(2)
plt.subplot(222)
plt.plot(fvec,db(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB)')

#het amplitudespectrum (lineair) als functie van het FFT-lijnnummer (DC is lijn 0)
plt.subplot(223)
plt.plot(np.abs(U),'*')
plt.xlabel('DFT lijnnummer')
plt.ylabel('Amplitude (lin)')

#het amplitudespectrum (lineair) als functie van de frequentie
plt.subplot(224)
plt.plot(fvec,np.abs(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (lin)')

# 3: sinus ? Hz; 
"""Herhaal deze oefening, maar kies nu een ! zodat een geheel aantal perioden gemeten wordt (bv.
één periode), hou fs en N constant. Wat merk je op?"""
    
A = 1
fs = 1000
N = 16
phi = np.pi/2

Tend = np.around(N*1/fs,5)
tvec = np.arange(0,Tend,1/fs)

f_NL = 1/Tend

u = np.sin(2*np.pi*f_NL*tvec+phi)    
U = np.fft.fft(u)/N

plt.figure(3)
plt.subplot(221)
plt.plot(tvec,u,'k*:')
plt.xlabel('Time (s)')
plt.ylabel('u(t)')

fvec = np.arange(0,fs,fs/N)

plt.figure(3)
plt.subplot(222)
plt.plot(fvec,db(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB)')


plt.subplot(223)
plt.plot(np.abs(U),'*')
plt.xlabel('DFT lijnnummer')
plt.ylabel('Amplitude (lin)')

plt.subplot(224)
plt.plot(fvec,np.abs(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (lin)')

# 4: sinus f_NL Hz N = 32
"""Herhaal de oefeningen nu (met de ‘goede’ !), maar verdubbel het aantal punten (N = 32),
maar houd Ts constant."""
    
A = 1
fs = 1000
N = 32
phi = np.pi/2

Tend = np.around(N*1/fs,5)
tvec = np.arange(0,Tend,1/fs)

u = np.sin(2*np.pi*f_NL*tvec+phi)    
U = np.fft.fft(u)/N

plt.figure(4)
plt.subplot(221)
plt.plot(tvec,u,'k*:')
plt.xlabel('Time (s)')
plt.ylabel('u(t)')

fvec = np.arange(0,fs,fs/N)

plt.figure(4)
plt.subplot(222)
plt.plot(fvec,db(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB)')


plt.subplot(223)
plt.plot(np.abs(U),'*')
plt.xlabel('DFT lijnnummer')
plt.ylabel('Amplitude (lin)')

plt.subplot(224)
plt.plot(fvec,np.abs(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (lin)')

# 5: sinus f_NL Hz 1 period N = 32
"""Je kan deze oefening nogmaals doen, waarbij je ten opzichte van oefening 3 opnieuw N verdubbelt
maar nu de eindtijd behoudt (1 periode). Interpreteer het verschil met hierboven."""

A = 1
N = 32
phi = np.pi/2

Tend = np.around(1/f_NL,5)
Ts = Tend/N
fs = 1/Ts
tvec = np.arange(0,Tend,Ts)

u = np.sin(2*np.pi*f_NL*tvec+phi)    
U = np.fft.fft(u)/N

plt.figure(5)
plt.subplot(221)
plt.plot(tvec,u,'k*:')
plt.xlabel('Time (s)')
plt.ylabel('u(t)')

fvec = np.arange(0,fs,fs/N)

plt.figure(5)
plt.subplot(222)
plt.plot(fvec,db(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB)')


plt.subplot(223)
plt.plot(np.abs(U),'*')
plt.xlabel('DFT lijnnummer')
plt.ylabel('Amplitude (lin)')

plt.subplot(224)
plt.plot(fvec,np.abs(U),'*')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (lin)')


# 6: Suppress leakage
"""Leakage is een representatie probleem van de data, geen probleem met de data aan zich. In
eerste instantie moet er steeds getracht worden een vensterbreedte te kiezen waarbinnen een
geheel aantal perioden valt. Lukt dat niet dan kan er overgegaan worden tot een ‘gewogen
venster’, t.t.z. een venster dat niet rechthoekig is. Het Hanning venster en de weging die het op
het signaal toepast zijn weergegeven in Fig. 3"""
fs = 1024
Ts = 1/fs
N = 1024

f1 = 20.5
f2 = 200.5

tvec = np.arange(0,np.around(N*Ts),Ts)
fvec = np.arange(0,fs,fs/N)

u = np.sin(2*np.pi*f1*tvec)+1e-4*np.sin(2*np.pi*tvec*f2)
U = np.fft.fft(u)/N

plt.figure(6)
plt.subplot(211)
plt.plot(fvec,db(U))
plt.xlabel('Frquency (Hz)')
plt.ylabel('|U| (dB)')
plt.xlim([0,fs/2])

window = np.hanning(N)
uHanning = np.multiply(u,window)

UHanning = np.fft.fft(uHanning)/N

plt.subplot(212)
plt.plot(fvec,db(UHanning))
plt.xlabel('Frquency (Hz)')
plt.ylabel('|U| (dB)')
plt.xlim([0,fs/2])

plt.figure(7)
plt.plot(tvec,uHanning)
plt.plot(tvec,u,':')
plt.plot(tvec,window,'r')
plt.xlabel('Time (s)')
plt.ylabel('u(n)')
plt.legend(('u(n)*window','u(n)','window'))
plt.show() 
