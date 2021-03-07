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
        self.Mass2 = 1.0
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
        Rsw[qi, :] = (1.0+inv*fsq(q, fs))**2*PolarW[qi, :]/denorm/inv

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

    return dRsT, dRaT


if __name__ == "__main__":

    Para = para()

    tgrid = np.loadtxt("tgrid.dat")
    kgrid = np.loadtxt("kgrid.dat")

    dRsT, dRaT = InterTau(tgrid, kgrid, Para)

    # print(polar.Polarisi(2.2e-16, 0.0, Para.EF))
    # print(dRsT[10, :])

    ########### Plot Polarization in Tau ################
    plt.figure()
    for qi, q in enumerate(kgrid[:10]):
        plt.plot(tgrid, dRsT[qi, :], label=f"{q}")
        # plt.plot(tgrid, dRsT[0, :], label=f"{q}")
    plt.title("dRs")
    plt.legend()
    plt.show()

    plt.figure()
    for qi, q in enumerate(kgrid[:10]):
        plt.plot(tgrid, dRaT[qi, :], label=f"{q}")
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
