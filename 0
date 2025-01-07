#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WPO2: Discrete-time signals, sampling and its issues. The theorem of Shannon. The DFT. (zie pdf opgave)

To learn:
    recognise aliasing in the time-domain
    recognise aliasing in the frequency-domain
    avoid aliasing

@author: ELV
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
# close all figures
plt.close('all')

# Oefening 1: Aliasing in the time-domain
endT = 1 #end time
f1 = 1  #freq
N_c = 1e3   #aantal pt vr discreet 
Ts_c = endT/N_c     #periode

f2 = 9
N_d = 10
Ts_d = endT/N_d

def mySin(f, tEnd, Ts): #maak sinus functie
    tVec = np.arange(0, tEnd,Ts)
    return(np.sin(2*np.pi*f*tVec), tVec)

y1c, t1c = mySin(1, 1, Ts_c)
y1d, t1d = mySin(1, 1, Ts_d)
y2c, t2c = mySin(9, 1, Ts_c)
y2d, t2d = mySin(9, 1, Ts_d)

plt.figure(1)
def plotter(sp, x1, y1, shape1, x2 = 1, y2 ="", shape2="" ):
    plt.subplot(sp)
    plt.plot(x1, y1, shape1)
    plt.plot(x2, y2, shape2)
    plt.xlabel('t(s)')
    plt.ylabel('y')
    
plotter(311, t1c, y1c, 'b', t1d, y1d, 'b*') #vergelijking sample rate on sinus1
plotter(312, t1d, y1d, 'b*', t2d, y2d, 'ro')    #vergelijking 2 discrete sinusen
plotter(313, t2c, y2c, 'r', t2d, y2d, 'ro') #vergelijking sample rate on sinus2
# sys.exit()
#%%
# Oefening 2: Aliasing looking for the Nyquist frequency
fs = 10 #bemonsteringsfrequentie in Hz, dit betekent dat de nyquist freq gelijk is aan fs/2 = 5
tvec = np.arange(0,1,1/fs) #ik ga dit niet gebruiken omdat ik in mySin de tijdsvector ook bepaal
N = len(tvec)
fs_c = 1e3
tc = np.arange(0,1,1/fs_c)


arrFreq = np.array([2,4,8,6]) #omdat 2Hz en 4Hz onder de Nyquist freq zitten gaat er hier geen aliasing voorkomen, we verwachten aliasing bij 8&6Hz

plt.figure(2)
for i in range(len(arrFreq)):
    y, t = mySin(arrFreq[i], 1, 1/fs)
    yc, tc = mySin(arrFreq[i], 1, 1/fs_c)

    plotter(220+i+1, t, y, 'r*', tc, yc, 'k')

#Herhaald voor een sinus van 5 Hz en 12 Hz   
y = np.sin(2*np.pi*5*tvec)
y_cont = np.sin(2*np.pi*5*tc)
plt.figure(3)
plt.title('5 Hz')
plt.plot(tvec,y,'r*:')
plt.plot(tc,y_cont,'k')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')

y = np.sin(2*np.pi*12*tvec)
y_cont = np.sin(2*np.pi*12*tc)

plt.figure(4)
plt.title('12 Hz')
plt.plot(tvec,y,'r*:')
plt.plot(tc,y_cont,'k')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
# sys.exit()
#%%  # Oefening 3: built-in DFT

"""Hergebruik de sinus met een frequentie van 3 Hz, gesampled aan 10 Hz, deze keer voor
een duurtijd van 2 s. Bereken hiervan het spectrum via een zelf geprogrammeerde
DFT-functie, en vergelijk dit dan met de ingebouwde functie. Denk eraan dat het
spectrum bestaat uit complexe getallen."""
fs = 10
y3, tvec = mySin(3, 2, 1/fs)

Y3 = np.fft.fft(y3)/1 #drawing DFT using np built in function
N = len(tvec)

# own DFT
Y3_own = np.zeros(Y3.shape[0],dtype=complex)
for k in range(N):
    for t in range(N):
        Y3_own[k] += y3[t]*np.exp(-2*1j*np.pi*k*t/N) #formule 2a, waarbij t-waardes n voorstellen

#plot amplitude in functie van k        
plt.figure(5)
plotter(111, np.arange(0,N), np.abs(Y3), 'ko-', np.arange(0,N) , np.abs(Y3_own), 'r*:')#vergelijk discretisatie eigen functie met built in functie
plt.xlabel('k')
plt.ylabel('|Y|')
plt.legend(('FFT','Own DFT'))       

# %% # own DFTon 2N lines
"""

Y3_own = np.zeros(2*Y3.shape[0],dtype=complex)
for k in range(2*N):
    for t in range(N):
        Y3_own[k] += y3[t]*np.exp(-2*1j*np.pi*k*t/(2*N))
        
plt.figure(8)
plt.plot(np.abs(Y3),'ko-')
plt.plot(np.abs(Y3_own),'r*:')
plt.xlabel('k')
plt.ylabel('|Y|')
plt.legend(('FFT','Own DFT'))"""
# %%

#plot amplitude in functie van frequentie
fvec = np.arange(0,fs,fs/N)

plt.figure(6)
plotter(111, fvec, np.abs(Y3),'ko-', fvec, np.abs(Y3_own),'r*:')
plt.xlabel('Frequency (Hz)')
plt.ylabel('|Y|')
plt.legend(('FFT','Own DFT'))
# sys.exit()

# %%  # Oefening 4

"""Herhaal nu oefeningen 2 uit de opdracht over Aliasing waarbij je nu het spectrum van de
verschillende sinussen vergelijkt.
Maak een vergelijking van sinussen met frequenties 2, 4, 6 en 8 Hz, nog steeds gebruik makend
van 10 punten en 1 s totale duur, en vergelijk met de ‘continue’ variant. Hoe manifesteert
Aliasing zich in het frequentiedomein?"""


fs = 10
tvec = np.arange(0,1,1/fs)
N = len(tvec)

fvec = np.arange(0,fs,fs/N)
#fvec = np.arange(-fs/2,fs/2,fs/N)

fs_cont = 1e3
tvec_cont = np.arange(0,1,1/fs_cont)


arrFreq = np.array([2,4,6,8])

plt.figure(7)
for i in range(len(arrFreq)):
    y = np.sin(2*np.pi*arrFreq[i]*tvec)
    y_cont = np.sin(2*np.pi*arrFreq[i]*tvec_cont)
    
    plt.subplot(4,2,(2*i)+1)
    plt.plot(tvec,y,'r*:')
    plt.plot(tvec_cont,y_cont,'k')
    plt.title(str(arrFreq[i]) + ' Hz')
#    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    
    plt.subplot(4,2,(2*i)+2)
    plt.stem(fvec,np.abs(np.fft.fft(y)),'r')
    plt.title(str(arrFreq[i]) + ' Hz')
#    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')

plt.subplot(4,2,7)
plt.xlabel('Time (s)') 

plt.subplot(4,2,8)
plt.xlabel('Frequency (Hz)')

"""conclusion: we notice that 4Hz & 8Hz,and 2Hz & 8Hz have the same spectral content, which
gives a wrong idea. hte higher frequenties need to be sampled at a higher rate""" 
