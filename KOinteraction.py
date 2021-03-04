#!/usr/bin/env python3
from utility.IO import *
from utility.plot import *
import utility.polar0 as polar
import utility.dlr_boson as dlr
import numpy as np
import scipy.linalg as slinalg
from scipy.linalg import lu_factor, lu_solve


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

    print("Maximum fitting error of Rs in frequency: ", np.max(abs(Rsw-_Rsw)))
    print("Maximum fitting error of Ra in frequency: ", np.max(abs(Raw-_Raw)))

    return dRsT, dRaT


if __name__ == "__main__":
    import grid
    Para = param("./")

    ############# Construct Tau and Mom grids ################
    T = grid.Tau()
    # double beta, int _size, double scale
    T.build(Para.Beta, Para.TauGridSize, 6.0/Para.EF)
    # print(T.grid)
    Kbose = grid.BoseK()
    # double kF, double maxK, int _size, double scale
    Kbose.build(Para.kF, Para.MaxExtMom,
                Para.MomGridSize, 1.0/Para.kF/2.0)
    # print(Kbose.grid)

    dRsT, dRaT = InterTau(T.grid, Kbose.grid, Para)

    ########### Plot Polarization in Tau ################
    plt.figure()
    for qi, q in enumerate(Kbose.grid[:10]):
        Errorbar(T.grid, dRsT[qi, :], label=f"{q}")
    plt.title("dRs")
    plt.legend()
    plt.show()

    plt.figure()
    for qi, q in enumerate(Kbose.grid[:10]):
        Errorbar(T.grid, dRaT[qi, :], label=f"{q}")
    plt.title("dRa")
    plt.legend()
    plt.show()

    with open("interaction.data", "w") as f:
        for qi, q in enumerate(Kbose.grid):
            for ti, t in enumerate(T.grid):
                f.write(f"{dRsT[qi, ti]} ")
        for qi, q in enumerate(Kbose.grid):
            for ti, t in enumerate(T.grid):
                f.write(f"{dRaT[qi, ti]} ")
