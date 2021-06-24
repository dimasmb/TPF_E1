# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 18:25:20 2021

@author: Dimas
"""
import numpy as np
import scipy.signal as ss
import matplotlib.pyplot as plt


def makeBode(num, den, z, p):

    if abs(z)<abs(p):
        if z==0:
            n=0
        else:
            n=int(np.log10(abs(z)))
        m=int(np.log10(abs(p)))
    else:
        if p==0:
            n=0
        else:
            n=int(np.log10(abs(p)))
        m=int(np.log10(abs(z)))

    H=ss.TransferFunction(num, den)
    w=np.logspace(n-3, m+3, 1000)
    bode=ss.bode(H, w=w)

    fig, (ax0, ax1)=plt.subplots(2)
    ax0.plot(w, bode[1])
    ax0.set_xscale("log")
    ax0.grid()

    ax0.set_title("Bode")

    if max(abs(bode[1]))<1:
        ax0.set_ylim(-1.2, 1.2)

    ax1.plot(w, bode[2])
    ax1.set_xscale("log")
    ax1.grid()
    ax1.set_title("Fase")
    if max(abs(bode[2]))<1:
        ax1.set_ylim(-1.2, 1.2)

    fig.tight_layout()


    return

def makeZP1(z, p):

    fig, ax=plt.subplots()
    ax.axhline(0, color='black')
    ax.axvline(0, color='black')
    ax.plot(z, 0, "o", label="ceros")
    ax.plot(p, 0, "x", label="polos", color='r')
    
    if abs(p)>abs(z):
        m=abs(p)
    else:
        m=abs(z)
        
    ax.set_xlim(-(m+m*0.1), m+m*0.1)
    plt.grid()
    plt.legend()

    
def pasalgo1(z, p, gain):
    
    num=[1*gain, -z*gain]
    den=[1, -p]

    makeBode(num, den, z, p)
    makeZP1(z, p)
    plt.show()
    return num, den


def pasalgo2(z, p, gain):

    num=[gain*1/(z[0]**2), gain*2*z[1]/z[0], 1*gain]
    den=[1/(p[0]**2), 2*p[1]/p[0], 1]

    makeBode(num, den, z[0], p[0])
    plt.show()
    return num, den


def pasaAltos(orden, z, p, gain):
    if orden==1:
        num=[gain, 0]
        den=[1/(-p), 1]
        z1=z
        p1=p
    else:
        num = [1*gain, 0, 0]
        den = [1 / (p[0] ** 2), 2 * p[1] / p[0], 1]
        z1 = z[0]
        p1 = p[0]
    makeBode(num, den, z1, p1)
    plt.show()
    return num, den


def pasaBajos(orden, z, p, gain):
    if orden==1:
        num=[gain]
        den=[1/(-p), 1]
        z1=z
        p1=p
    else:
        num = [gain*1]
        den = [1 / (p[0] ** 2), 2 * p[1] / p[0], 1]
        z1 = z[0]
        p1 = p[0]
    makeBode(num, den, z1, p1)
    plt.show()
    return num, den

def pasaTodo(orden, z, p, gain):
    if orden==1:
        num=[1/(-p), -1]
        den=[1/(-p), 1]
        z1 = z
        p1 = p
    else:
        num = [gain * 1 / (z[0] ** 2), -gain * 2 * z[1] / z[0], 1 * gain]
        den = [1 / (p[0] ** 2), 2 * p[1] / p[0], 1]
        z1 = z[0]
        p1 = p[0]
    makeBode(num, den, z1, p1)
    plt.show()
    return num, den

def pasaBanda(z, p, gain):
    num = [0, gain * 1, 0]
    den = [1 / (p[0] ** 2), 2 * p[1] / p[0], 1]
    makeBode(num, den, z[0], p[0])
    plt.show()
    return num, den

def notch(z, p, gain):
    num = [gain * 1 / (z[0] ** 2), 0, 1 * gain]
    den = [1 / (p[0] ** 2), 2 * p[1] / p[0], 1]
    makeBode(num, den, z[0], p[0])
    plt.show()
    return num, den

def RLC(R, L, C, Vo, nodos):
    if C==0:
        if Vo=="R":
            num=[R]
        elif Vo=="L":
            num=[L, 0]
        elif Vo == "C":
            num = [1]
        den=[L, R]
        corte=int(np.log10(abs(R/L)))
    else:
        if Vo=="R":
            num=[R*C, 0]
        elif Vo=="L":
            num=[C*L, 0, 0]
        elif Vo == "C":
            num = [1]
        den=[L*C, R*C, 1]
        corte=int(np.log10(np.sqrt(1/(abs(L*C)))))


    H = ss.TransferFunction(num, den)
    w = np.logspace(corte - 3, corte + 3, 1000)
    bode = ss.bode(H, w=w)

    fig, (ax0, ax1) = plt.subplots(2)
    ax0.plot(w, bode[1])
    ax0.set_xscale("log")
    ax0.grid()
    tit="Bode RLC " + nodos
    ax0.set_title(tit)

    ax1.plot(w, bode[2])
    ax1.set_xscale("log")
    ax1.grid()
    ax1.set_title("Fase")
    fig.tight_layout()
    plt.show()
    return


def parseSci(str):
    mult=str[-1]

    if mult.isnumeric():
        try:
            val = float(str)
        except ValueError:
            return "E"
    elif  mult=='p' or mult=='n' or mult=='u' or mult=='m' or mult=='k' or mult=='M' or mult=='G':
        try:
            val=float(str[:-1])
        except ValueError:
            return "E"
    else:
        return "E"


    if mult=='p':
        val=val*(10**(-12))
    elif mult=='n':
        val=val*(10**(-9))
    elif mult=='u':
        val=val*(10**(-6))
    elif mult=='m':
        val=val*(10**(-3))
    elif mult=='k':
        val=val*(10**(3))
    elif mult=='M':
        val=val*(10**(6))
    elif mult=='G':
        val=val*(10**(9))

    return val

def parseArray(arr):
    if arr.find(",") !=-1:
        strs=arr.split(",")
        val1=parseSci(strs[0])
        val2 = parseSci(strs[1])
        if val1=="E" or val2=="E":
            return "E"

        return val1, val2
    else:
        val=parseSci(arr)
        if val=="E":
            return "E"
        return val