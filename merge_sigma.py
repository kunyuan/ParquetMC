import fourier
from utility import *
from grid import *


def SigmaStatic(data, norm, Para):
    order1 = [np.average(d[1, :, :], axis=1) for d in data]
    static, staticErr = Estimate(order1, norm)
    return static, staticErr


def SigmaT(data, norm, Para):
    sumOrders = [np.sum(d[2:Para.Order+1, ...], axis=0) for d in data]
    dynamic, dynErr = Estimate(sumOrders, norm)
    return dynamic, dynErr


def PlotSigmaW(SigmaT, idx, Save=True):
    SigmaW, Spec = Fourier.SpectralT2W(SigmaT[idx, :])
    C0 = SigmaT[idx, 0]+SigmaT[idx, -1]
    SigmaWp = Fourier.naiveT2W(SigmaT[idx, :])

    fig, (ax1, ax2) = plt.subplots(1, 2)
    # ax.plot(Freq, Spec[idx, :], "k--",
    #         label="k={0}, spectral".format(Grid.MomGrid[idx]))
    k = Grid.MomGrid[idx]/Para.kF
    ax1.plot(phyFreq, SigmaW.imag, "ro",
             label="$k/k_F={0}$, spectral".format(k))
    ax1.plot(phyFreq, SigmaWp.imag, "b--",
             label="$k/k_F={0}$, naive".format(k))
    ax1.plot(phyFreq, -C0/phyFreq, "g--",
             label="$k/k_F={0}$, tail".format(k))

    # ax1.plot(phyFreq, SigmaW.real, "ro",
    #          label="$k/k_F={0}$, spectral".format(k))
    # ax1.plot(phyFreq, SigmaWp.real, "b--",
    #          label="$k/k_F={0}$, naive".format(k))

    ax1.set_ylim([-0.12, 0.12])

    SigmaTp = Fourier.naiveW2T(SigmaW)
    ax2.plot(Grid.TauGrid, Dynamic[idx, :], "r--",
             label="$k/k_F={0}$, original".format(k))
    ax2.plot(Grid.TauGrid, SigmaTp.real, "b--",
             label="$k/k_F={0}$, spectral fourier".format(k))

    ax1.legend(loc=1, frameon=False, fontsize=size)
    ax2.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("SigmaW.pdf")
    else:
        plt.show()


def PlotDataK(dataK, wn, FreqGrid, Save=True):
    fig, ax = plt.subplots()
    ax.plot(Grid.MomGrid, dataK[:, wn].real, "ro-",
            label="$n={0}$, real".format(FreqGrid[wn]))
    ax.plot(Grid.MomGrid, dataK[:, wn].imag, "bo-",
            label="$n={0}$, imag".format(FreqGrid[wn]))

    plt.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("DataK.pdf")
    else:
        plt.show()


def PlotSigmaT(dataW, idx, Save=True):
    dataT, Spec = Fourier.SpectralW2T(dataW[idx, :])
    dataTp = Fourier.naiveW2T(dataW[idx, :])

    k = Grid.MomGrid[idx]/Para.kF
    fig, ax = plt.subplots()
    ax.plot(Grid.TauGrid, dataT.real, "r--",
            label="$k/k_F={0}$, spectral".format(k))
    ax.plot(Grid.TauGrid, dataTp.real, "b--",
            label="$k/k_F={0}$, naive".format(k))

    # ax1.set_ylim([-0.12, 0.12])

    plt.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("DeltaG.pdf")
    else:
        plt.show()


if __name__ == "__main__":

    Para = param()
    Grid = grid(Para)
    Spectral = spectral(Para)
    Order = range(0, Para.Order+1)

    MaxFreq = 1024
    Freq = np.array(range(-MaxFreq, MaxFreq))
    # # print len(Freq)
    phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency

    Fourier = fourier.fourier(Grid.TauGrid, phyFreq, Para.Beta)
    Fourier.InitializeKernel(100.0, 1024, "Fermi", 1.0e-13)

    folder = "./Data"

    filename = "sigma_pid[0-9]+.dat"

    shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)

    Data, Norm, Step = LoadFile(folder, filename, shape)

    Static, StaticErr = SigmaStatic(Data, Norm, Para)
    Dynamic, DynErr = SigmaT(Data, Norm, Para)

    print "Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr))

    # PlotSigmaW(Dynamic, int(Para.MomGridSize/3), True)

    SigmaW, _ = Fourier.SpectralT2W(Dynamic)

    PlotDataK(SigmaW, MaxFreq, Freq, False)

    BareG = np.zeros((Grid.MomSize, len(phyFreq)), dtype=complex)
    for i, q in enumerate(Grid.MomGrid):
        BareG[i, :] = -1.0/(1j*phyFreq-q*q+Para.EF)

    dG_W = 1.0/(1.0/BareG-SigmaW)
    dG_W = SigmaW*BareG*BareG/(1-SigmaW*BareG)

    dG_T, _ = Fourier.SpectralW2T(dG_W)
    dG_Tp = Fourier.naiveW2T(dG_W)

    print "Maximum Error of \delta G: ", np.amax(abs(dG_T-dG_Tp))

    PlotSigmaT(dG_W, int(Para.MomGridSize/3), True)
