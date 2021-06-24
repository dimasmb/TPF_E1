import numpy as np
import scipy.signal as ss
import matplotlib.pyplot as plt
import resources.pasa as pa


def salida(tipo, amp, num, den):
    if tipo=="Escalon":

        # t, y=ss.step((round(amp)*num, den), T=t)
        n=np.array(num)
        d=np.array(den)

        n=n*amp
        d = np.append(d, 0)
        H = ss.TransferFunction(n, d)
        w = np.logspace(-3, 7, 1000)
        # w = np.logspace(n - 3, m + 3, 1000)
        y = ss.bode(H, w=w)


        fig, (ax0, ax1) = plt.subplots(2)
        ax0.plot(w, y[1])
        ax0.set_xscale("log")
        ax0.grid()
        ax0.set_title("Bode")

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
    elif tipo=="Seno":

        pass
    return