#!/usr/bin/env python3
from utility.IO import *
from utility.plot import *
import utility.polar0 as polar
import utility.dlr_boson as dlr
import numpy as np
import scipy.linalg as slinalg
from scipy.linalg import lu_factor, lu_solve
import utility.angle as legendre
from scipy import integrate


def fsq(q, fs): return fs
def faq(q, fa): return fa


def InterFreq(wngrid, kgrid, Para):
    """wngrid: matsubara frequency (integers!)
        dRsW=(Physical dRsW)/v_q
    """
    ksize = len(kgrid)
    wnsize = len(wngrid)
    PolarW = np.zeros([ksize, wnsize])
    wngrid = np.array(wngrid)*2.0*np.pi/Para.Beta
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
        Rsw[qi, :] = (1.0+inv*fsq(q, fs))**2*PolarW[qi, :]/denorm

        # dRa=fa^2*Pi0/(1-fa*Pi0)
        denorm = 1.0-faq(q, fa)*PolarW[qi, :]
        Raw[qi, :] = faq(q, fa)**2*PolarW[qi, :]/denorm

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

    Rsw, Raw = InterFreq(wnGrid, kgrid, Para)
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

def LandauParameter(Para):
    theta = np.linspace(0.0, np.pi, 1000)
    Qgrid = Para.kF*2.0*np.sin(theta/2.0)

    dRsW, dRaW = InterFreq([0, ], Qgrid, Para)

    dRsW = dRsW[:, 0]
    dRsW = (1.0+dRsW)*8.0*np.pi/(Qgrid**2+Para.Mass2)
    Rs0 = legendre.LegendreCoeff(dRsW, -np.cos(theta), [0, ], 0)[0]
    # print(dRsW*Para.Nf)
    # print("l=0: ", Rs0*Para.Nf)
    Fs=-0.5*Rs0*Para.Nf
    Fa=-0.5*Rs0*Para.Nf
    return Fs, Fa



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

    # tgrid = np.loadtxt("./calculator/tgrid.dat")
    # kgrid = np.loadtxt("./calculator/kgrid.dat")
    # dRsT, dRaT = InterTau(tgrid, kgrid, Para)

    dRsT, dRaT = InterTau(T.grid, Kbose.grid, Para)

    Fs, Fa=LandauParameter(Para)
    print(Fs, Fa)


    ########### Plot Polarization in Tau ################
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    N = 5
    for i in range(N):
        qi = int(i*len(Kbose.grid)/N)
        ax1.plot(T.grid, dRsT[qi, :], label=f"{Kbose.grid[qi]}")
    # for qi, q in enumerate(Kbose.grid[::20]):
    #     Errorbar(T.grid, dRsT[qi, :], label=f"{q}")
    ax1.set_title("dRs")

    ax2.plot(np.array(Kbose.grid)/Para.kF,
             dRsT[:, 0], label=f"$\\tau={T.grid[0]}$")
    ax2.set_xlabel("q/kF")
    ax2.set_title("dRs in K")

    N = 5
    for i in range(N):
        qi = int(i*len(Kbose.grid)/N)
        ax3.plot(T.grid, dRaT[qi, :], label=f"{Kbose.grid[qi]}")
    # for qi, q in enumerate(Kbose.grid[::20]):
    #     Errorbar(T.grid, dRsT[qi, :], label=f"{q}")
    ax3.set_title("dRs")

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()

    with open("interaction.data", "w") as f:
        for qi, q in enumerate(Kbose.grid):
            for ti, t in enumerate(T.grid):
                f.write(f"{dRsT[qi, ti]} ")
        for qi, q in enumerate(Kbose.grid):
            for ti, t in enumerate(T.grid):
                f.write(f"{dRaT[qi, ti]} ")
