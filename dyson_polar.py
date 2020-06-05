#!/usr/bin/env python3
from utility.IO import *
from utility.plot import *
import utility.fourier as fourier


def PlotPolarW(PolarT, MomGrid, idx, Save=True):
    PolarW, Spec = Fourier.SpectralT2W(PolarT[idx, :])
    # C0 = SigmaT[idx, 0]+SigmaT[idx, -1]
    PolarWp = Fourier.naiveT2W(PolarT[idx, :])

    _, (ax1, ax2) = plt.subplots(1, 2)
    # ax.plot(Freq, Spec[idx, :], "k--",
    #         label="k={0}, spectral".format(Grid.MomGrid[idx]))
    k = MomGrid[idx]/Para.kF
    ax1.plot(phyFreq, PolarW.real, "ro", markersize=3,
             label="$k/k_F={0:.4f}$, spectral".format(k))
    ax1.plot(phyFreq, PolarWp.real, "b--",
             label="$k/k_F={0:.4f}$, naive".format(k))
    # ax1.plot(phyFreq, -C0/phyFreq, "g--",
    #          label="$k/k_F={0:.4f}$, tail".format(k))
    # ax1.plot(phyFreq, (-1/(1.0j*phyFreq)).imag, "g--",
    #          label="$k/k_F={0:.4f}$, tail".format(k))

    # ax1.plot(phyFreq, PolarW.real, "bo",
    #          label="$k/k_F={0}$, spectral".format(k))
    # ax1.plot(phyFreq, PolarWp.real, "b--",
    #          label="$k/k_F={0}$, naive".format(k))

    ax1.set_ylim([-0.12, 0.12])

    PolarTspectral, _ = Fourier.SpectralW2T(PolarW)
    PolarTnaive = Fourier.naiveW2T(PolarWp)
    ax2.plot(TauGrid, Dynamic[idx, :], "k--",
             label="$k/k_F={0:.4f}$, original".format(k))
    ax2.plot(TauGrid, PolarTspectral.real, "r--",
             label="$k/k_F={0:.4f}$, spectral fourier".format(k))
    ax2.plot(TauGrid, PolarTnaive.real, "b--",
             label="$k/k_F={0:.4f}$, naive fourier".format(k))

    ax1.legend(loc=1, frameon=False, fontsize=size)
    ax2.legend(loc=1, frameon=False, fontsize=size)
    plt.grid()
    if Save:
        plt.savefig("PolarW.pdf")
    else:
        plt.show()


if __name__ == "__main__":

    Para = param()
    Order = range(0, Para.Order+1)

    MaxFreq = 1024
    Freq = np.array(range(0, MaxFreq))
    phyFreq = (Freq*2.0)*np.pi/Para.Beta  # the physical frequency

    shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)
    Data, Norm, Step, Grids = LoadFile("./Data", "polar_pid[0-9]+.dat", shape)

    TauGrid = Grids["TauGrid"]
    MomGrid = Grids["KGrid"]

    Fourier = fourier.fourierBose(TauGrid, phyFreq, Para.Beta)
    Fourier.InitializeKernel(100.0, 1024, "Bose", 1.0e-13)

    # Dynamic, DynErr = Estimate(
    #     Data, Norm, lambda d: np.sum(d[1:Para.Order+1, ...], axis=0))
    Dynamic, DynErr = Estimate(
        Data, Norm, lambda d: np.sum(d[1:2, ...], axis=0))

    PlotPolarW(Dynamic, MomGrid, 0, False)
    # PlotDataK(SigmaW, MaxFreq, Freq, False)
