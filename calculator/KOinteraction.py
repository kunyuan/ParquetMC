#!/usr/bin/env python3
import dlr_boson as dlr
import polar0 as polar
import numpy as np
import scipy.linalg as slinalg
from scipy.linalg import lu_factor, lu_solve
import matplotlib.pyplot as plt


class para:
    def __init__(self):
        self.Beta = 25.0
        self.Rs = 1.0
        self.Mass2 = 0.01
        self.Fs = 0.0
        self.Fa = 0.0
        self.kF = (9.0*np.pi/4.0)**(1.0/3.0)/self.Rs
        self.EF = self.kF**2
        self.Spin = 2
        self.Nf = self.kF/4.0/np.pi**2*self.Spin
        self.Beta = self.Beta/self.EF


def fsq(q, fs): return fs
def faq(q, fa): return fa


def InterFreq(wngrid, kgrid, Para, addBare):
    """wngrid: matsubara frequency (integers!)
    """
    ksize = len(kgrid)
    wnsize = len(wngrid)
    PolarW = np.zeros([ksize, wnsize])
    wngrid = wngrid*2.0*np.pi/Para.Beta
    for qi, q in enumerate(kgrid):
        for wi, w in enumerate(wngrid):
            PolarW[qi, wi] = polar.Polarisi(q, w, Para.EF)
            print(q, w, PolarW[qi, wi])

    # print(PolarW[:, 0])

    Rsw = np.zeros_like(PolarW)
    Raw = np.zeros_like(PolarW)

    fs = Para.Fs/Para.Nf
    fa = Para.Fa/Para.Nf

    for qi, q in enumerate(kgrid):
        # add a small mass to regularize everything
        inv = (q*q+Para.Mass2)/(8.0*np.pi)
        # dRs=(v+fs)^2*Pi0/(1-(v+fs)*Pi0)
        denorm = inv-(1.0+inv*fsq(q, fs))*PolarW[qi, :]
        # Rsw[qi, :] = (1.0+inv*fsq(q, fs))**2*PolarW[qi, :]/denorm/inv
        Rsw[qi, :] = (1.0+inv*fsq(q, fs))**2*PolarW[qi, :]/denorm

        # dRa=fa^2*Pi0/(1-fa*Pi0)
        denorm = 1.0-faq(q, fa)*PolarW[qi, :]
        Raw[qi, :] = faq(q, fa)**2*PolarW[qi, :]/denorm

        if addBare:
            Rsw[qi, :] += 1.0/inv

    return Rsw, Raw


def InterTau(tgrid, kgrid, Para, eps=1.0e-12):
    """only include the tau dependent part of the KO interaction"""

    ############  DLR basis   #################################
    Wmax = 10.0
    Nw = 1024
    Nt = int(Para.Beta*Wmax*100)
    Nwn = int(Wmax*Para.Beta/2.0/np.pi*100)
    # print("Nt: ", Nt)
    # print("Nwn: ", Nwn)
    k, wGrid, wnGrid, _tGrid = dlr.getDLR(Para.Beta, Wmax, Nw, Nwn, Nt, eps)

    KerW = dlr.getKerW(wnGrid, wGrid, Para.Beta)
    lu, piv = lu_factor(KerW)

    KerT = dlr.getKerT(tgrid, wGrid, Para.Beta)

    Rsw, Raw = InterFreq(wnGrid, kgrid, Para, False)
    print(Rsw[0, 0]+8.0*np.pi/Para.Mass2-8.0 *
          np.pi/(Para.Mass2+8.0*np.pi*Para.Nf))

    _Rsw, _Raw = np.zeros_like(Rsw), np.zeros_like(Raw)
    dRsT = np.zeros([len(kgrid), len(tgrid)])
    dRaT = np.zeros_like(dRsT)

    coeff = np.zeros([len(kgrid), k])

    for qi, _ in enumerate(kgrid):
        # dRs=(v+fs)^2*Pi0/(1-(v+fs)*Pi0)
        coeff[qi, :] = lu_solve((lu, piv), Rsw[qi, :])
        _Rsw[qi, :] = KerW @ coeff[qi, :]
        dRsT[qi, :] = KerT @ coeff[qi, :]

        # dRa=fa^2*Pi0/(1-fa*Pi0)
        coeff[qi, :] = lu_solve((lu, piv), Raw[qi, :])
        _Raw[qi, :] = KerW @ coeff[qi, :]
        dRaT[qi, :] = KerT @ coeff[qi, :]

    # print(np.average(dRsT[0, :])*Para.Beta, Rsw[0, 0])
    print("RsT(q=0, tau)-Rs(q=0, w=0)/beta: ",
          np.max(abs(dRsT[0, :]-Rsw[0, 0]/Para.Beta)))
    # print(dRsT[0, :])

    print("Maximum fitting error of Rs in frequency: ", np.max(abs(Rsw-_Rsw)))
    print("Maximum fitting error of Ra in frequency: ", np.max(abs(Raw-_Raw)))

    # q = 1.0*Para.kF
    # Nk = int(len(kgrid)/2)
    # q = kgrid[Nk]
    # wt = np.zeros(len(tgrid))
    # wn = np.array(range(-100, 101))*2*np.pi/Para.Beta

    # p = [polar.Polarisi(q, abs(w), Para.EF) for w in wn]
    # p = np.array(p)
    # w = 8.0*np.pi*p/(q*q+Para.Mass2-8.0*np.pi*p)*8.0*np.pi/(q*q+Para.Mass2)
    # print(w)
    # print(np.sum(w)/Para.Beta)
    # print(polar.Polarisi(1.0e-4, 0.0, Para.EF))
    # print(8.0*np.pi*0.097/(0.1+8.0*np.pi*0.097)
    #       * 8.0*np.pi/(q*q+Para.Mass2)/Para.Beta)
    # wt[ti] = np.sum(w*np.exp(1j*wn*t))/Para.Beta

    # for ti, t in enumerate(tgrid):
    #     p = [polar.Polarisi(q, abs(w), Para.EF) for w in wn]
    #     p = np.array(p)
    #     v = 8.0*np.pi/(q*q+Para.Mass2)
    #     Wint = 8.0*np.pi*p/(q*q+Para.Mass2-8.0*np.pi*p) * v
    #     wt[ti] = np.sum(Wint*np.cos(wn*t))/Para.Beta
    # print(p)
    # print(Wint)
    # print(Wint*np.cos(wn*t))
    # print(t, wt[ti])

    # print(wt[0])
    # plt.figure()
    # plt.plot(tgrid, wt)
    # plt.plot(tgrid, dRsT[Nk, :])
    # plt.show()

    return dRsT, dRaT


def RinT(tgrid, q):
    wn = np.array(range(-100, 101))*2*np.pi/Para.Beta
    p = np.array([polar.Polarisi(q, abs(w), Para.EF) for w in wn])
    v = 8.0*np.pi/(q*q+Para.Mass2)
    rint = np.zeros(len(tgrid))
    Wint = 8.0*np.pi*p/(q*q+Para.Mass2-8.0*np.pi*p) * v

    for ti, t in enumerate(tgrid):
        rint[ti] = np.sum(Wint*np.cos(wn*t))/Para.Beta

    return rint


if __name__ == "__main__":

    Para = para()

    tgrid = np.loadtxt("tgrid.dat")
    kgrid = np.loadtxt("kgrid.dat")

    dRsT, dRaT = InterTau(tgrid, kgrid, Para)

    # print(polar.Polarisi(2.2e-16, 0.0, Para.EF))
    # print(dRsT[10, :])
    print(kgrid)
    print(kgrid[-10::2])

    ########### Plot Polarization in Tau ################
    plt.figure()

    N = 5
    # for i in range(N):
    #     qi = int(i*len(kgrid)/N)
    #     plt.plot(tgrid, dRsT[qi, :], label=f"{kgrid[qi]}")
    #     plt.plot(tgrid, RinT(tgrid, kgrid[qi]), label=f"{kgrid[qi]}")

    qi = 0
    plt.plot(tgrid, dRsT[qi, :], label=f"{kgrid[qi]}")
    plt.plot(tgrid, RinT(tgrid, kgrid[qi]), label=f"{kgrid[qi]}")

    qi = -180
    plt.plot(tgrid, dRsT[qi, :], label=f"{kgrid[qi]}")
    plt.plot(tgrid, RinT(tgrid, kgrid[qi]), label=f"{kgrid[qi]}")

    plt.title("dRs")
    plt.legend()
    plt.show()

    plt.plot(kgrid/Para.kF, dRsT[:, 0], label=f"{tgrid[0]}")
    plt.xlabel("q/kF")
    plt.title("dRs in K")
    plt.legend()
    plt.show()

    plt.figure()
    for i in range(N):
        qi = int(i*len(kgrid)/N)
        plt.plot(tgrid, dRaT[qi, :], label=f"{kgrid[qi]}")
        # plt.plot(tgrid, RinT(tgrid, q), label=f"{q}")
    plt.title("dRa")
    plt.legend()
    plt.show()

    # savetxt("interaction.dat")
    np.save("Rs", dRsT)
    np.save("Ra", dRaT)

    # with open("interaction.dat", "w") as f:
    #     for qi, q in enumerate(kgrid):
    #         for ti, t in enumerate(tgrid):
    #             f.write(f"{dRsT[qi, ti]} ")
    # f.write("\n")
    # for qi, q in enumerate(kgrid):
    #     for ti, t in enumerate(tgrid):
    #         f.write(f"{dRaT[qi, ti]} ")
