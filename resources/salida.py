import numpy as np
import scipy.signal as ss
import matplotlib.pyplot as plt
import resources.pasa as pa


def salida(tipo, amp, num, den, frec=0):
    n = np.array(num)
    d = np.array(den)
    if tipo=="Escalon":
        n=n*amp
        d = np.append(d, 0)
        fc=int(np.sqrt(abs(1/d[0])))
        mini=int(np.log10(fc))
        maxi=mini
    else:   #if tipo=="Seno":
        ord=len(den)
        if type(amp)!=list:
            return "E"
        ampli=amp[0]
        frec=amp[1]
        a=2*np.pi*frec
        n=n*amp
        if ord=="1":
            d=np.append(d, [(a**2)*n[0], a])
        else:
            d=np.array([d[0], d[1], d[0]*(a**2)+1, d[1]*(a**2), a**2])
        fc = int(np.sqrt(abs(1 / d[0])))
        mini=min(fc, a)
        maxi=max(fc, a)


    H = ss.TransferFunction(n, d)
    w = np.logspace(mini - 3, maxi + 3, 1000)
    y = ss.bode(H, w=w)

    fig, (ax0, ax1) = plt.subplots(2)
    ax0.plot(w, y[1])
    ax0.set_xscale("log")
    ax0.grid()
    st="Bode salida con entrada " + tipo.lower()
    ax0.set_title(st)

    if max(abs(y[1])) < 1:
        ax0.set_ylim(-1.2, 1.2)

    ax1.plot(w, y[2])
    ax1.set_xscale("log")
    ax1.grid()
    ax1.set_title("Fase")
    if max(abs(y[2])) < 1:
        ax1.set_ylim(-1.2, 1.2)

    fig.tight_layout()
    plt.show()
    return