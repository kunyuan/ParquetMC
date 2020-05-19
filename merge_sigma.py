from utility.IO import *
import utility.fourier as fourier


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

    _, (ax1, ax2) = plt.subplots(1, 2)
    # ax.plot(Freq, Spec[idx, :], "k--",
    #         label="k={0}, spectral".format(Grid.MomGrid[idx]))
    k = MomGrid[idx]/Para.kF
    ax1.plot(phyFreq, SigmaW.imag, "ro", markersize=3,
             label="$k/k_F={0:.4f}$, spectral".format(k))
    ax1.plot(phyFreq, SigmaWp.imag, "b--",
             label="$k/k_F={0:.4f}$, naive".format(k))
    ax1.plot(phyFreq, -C0/phyFreq, "g--",
             label="$k/k_F={0:.4f}$, tail".format(k))

    # ax1.plot(phyFreq, SigmaW.real, "ro",
    #          label="$k/k_F={0}$, spectral".format(k))
    # ax1.plot(phyFreq, SigmaWp.real, "b--",
    #          label="$k/k_F={0}$, naive".format(k))

    ax1.set_ylim([-0.12, 0.12])

    SigmaTp = Fourier.naiveW2T(SigmaW)
    ax2.plot(TauGrid, Dynamic[idx, :], "r--",
             label="$k/k_F={0:.4f}$, original".format(k))
    ax2.plot(TauGrid, SigmaTp.real, "b--",
             label="$k/k_F={0:.4f}$, spectral fourier".format(k))

    ax1.legend(loc=1, frameon=False, fontsize=size)
    ax2.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("SigmaW.pdf")
    else:
        plt.show()


def PlotDataK(dataK, wn, FreqGrid, Save=True):
    fig, ax = plt.subplots()
    ax.plot(MomGrid, dataK[:, wn].real, "ro-",
            label="$n={0}$, real".format(FreqGrid[wn]))
    ax.plot(MomGrid, dataK[:, wn].imag, "bo-",
            label="$n={0}$, imag".format(FreqGrid[wn]))

    plt.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("DataK.pdf")
    else:
        plt.show()


def PlotSigmaT(dataW, kList, Save=True):
    fig, ax = plt.subplots()

    for idx in kList:
        k = MomGrid[idx]/Para.kF
        dataT, _ = Fourier.SpectralW2T(dataW[idx, :])
        dataTp = Fourier.naiveW2T(dataW[idx, :])
        ax.plot(TauGrid, dataT.real, "--",
                label="$k/k_F={0:.4f}$, spectral".format(k))
        # ax.plot(TauGrid, dataTp.real, "--",
        #         label="$k/k_F={0}$, naive".format(k))

    # ax1.set_ylim([-0.12, 0.12])

    plt.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("DeltaG.pdf")
    else:
        plt.show()


if __name__ == "__main__":

    Para = param()
    Order = range(0, Para.Order+1)

    MaxFreq = 1024
    Freq = np.array(range(-MaxFreq, MaxFreq))
    phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency

    shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)
    Data, Norm, Step, Grids = LoadFile("./Data", "sigma_pid[0-9]+.dat", shape)

    TauGrid = Grids["TauGrid"]
    MomGrid = Grids["KGrid"]

    Fourier = fourier.fourier(TauGrid, phyFreq, Para.Beta)
    Fourier.InitializeKernel(100.0, 1024, "Fermi", 1.0e-13)

    Static, StaticErr = SigmaStatic(Data, Norm, Para)
    Dynamic, DynErr = SigmaT(Data, Norm, Para)

    print "Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr))

    PlotSigmaW(Dynamic, int(Para.MomGridSize/3), False)

    SigmaW, _ = Fourier.SpectralT2W(Dynamic)

    PlotDataK(SigmaW, MaxFreq, Freq, False)

    BareG = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)
    for i, q in enumerate(MomGrid):
        BareG[i, :] = -1.0/(1j*phyFreq-q*q+Para.EF)

    dG_W = 1.0/(1.0/BareG-SigmaW)
    dG_W = SigmaW*BareG*BareG/(1-SigmaW*BareG)

    dG_T, _ = Fourier.SpectralW2T(dG_W)
    dG_Tp = Fourier.naiveW2T(dG_W)

    print "Maximum Error of \delta G: ", np.amax(abs(dG_T-dG_Tp))

    PlotSigmaT(dG_W, range(0, Para.MomGridSize, Para.MomGridSize/8), False)

    with open("green.dat", "w") as f:
        for k in range(Para.MomGridSize):
            f.write("{0} ".format(Static[k]))
        f.write("\n")

        for k in range(Para.MomGridSize):
            for t in range(Para.TauGridSize):
                f.write("{0} ".format(dG_T[k, t].real))
        f.write("\n")
