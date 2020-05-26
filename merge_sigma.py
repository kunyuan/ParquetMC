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


def PlotSigmaW(SigmaT, MomGrid, idx, Save=True):
    SigmaW, Spec = Fourier.SpectralT2W(SigmaT[idx, :])
    C0 = SigmaT[idx, 0]+SigmaT[idx, -1]
    SigmaWp = Fourier.naiveT2W(SigmaT[idx, :])

    _, (ax1, ax2) = plt.subplots(1, 2)
    # ax.plot(Freq, Spec[idx, :], "k--",
    #         label="k={0}, spectral".format(Grid.MomGrid[idx]))
    k = MomGrid[idx]/Para.kF
    ax1.plot(phyFreq, SigmaW.imag, "ro", markersize=3,
             label="$k/k_F={0:.4f}$, spectral".format(k))
    # ax1.plot(phyFreq, SigmaWp.imag, "b--",
    #          label="$k/k_F={0:.4f}$, naive".format(k))
    ax1.plot(phyFreq, -C0/phyFreq, "g--",
             label="$k/k_F={0:.4f}$, tail".format(k))
    # ax1.plot(phyFreq, (-1/(1.0j*phyFreq)).imag, "g--",
    #          label="$k/k_F={0:.4f}$, tail".format(k))

    ax1.plot(phyFreq, SigmaW.real, "bo",
             label="$k/k_F={0}$, spectral".format(k))
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

    plt.axvline(x=Para.kF, linestyle='--')

    plt.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("DataK.pdf")
    else:
        plt.show()


def PlotSigmaT(dataW, kList, Save=True):
    fig, ax = plt.subplots()

    for idx in kList:
        k = MomGrid[idx]
        dataT, _ = Fourier.SpectralW2T(dataW[idx, :])
        dataTp = Fourier.naiveW2T(dataW[idx, :])
        ax.plot(TauGrid/Para.Beta, dataT.real, "--",
                label="$k/k_F={0:.4f}$, spectral".format(k/Para.kF))
        # ax.plot(TauGrid, dataTp.real, "--",
        #         label="$k/k_F={0}$, naive".format(k))
        Ek = k*k-Para.EF
        # print "Ek: ", Ek
        # G0 = np.exp(-Ek*TauGrid)/(1.0+np.exp(-Ek*Para.Beta))
        # ax.plot(TauGrid/Para.Beta, G0, "--", label="Bare G")

    ax.set_xlim([TauGrid[0]/Para.Beta, TauGrid[-1]/Para.Beta])

    plt.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("DeltaG.pdf")
    else:
        plt.show()


def PlotG(dataT, kFidx, Save=True):
    _, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(MomGrid/Para.kF, dataT[:, -1], "ro-",
             label="$\\tau=-0^+$")
    ax1.plot(MomGrid/Para.kF, dataT[:, 0], "bs-",
             label="$\\tau=0^+$")
    # ax.plot(TauGrid, dataTp.real, "--",
    #         label="$k/k_F={0}$, naive".format(k))

    ax1.plot([MomGrid[0], MomGrid[-1]], [0.5, 0.5], linestyle='--')
    ax1.axvline(x=1.0, linestyle='--')
    ax1.set_xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    ax1.set_xlabel("$k/k_F$")
    # ax1.set_ylim([-0.12, 0.12])

    for idx in [kFidx/2, kFidx, kFidx*2]:
        k = MomGrid[idx]
        ax2.plot(TauGrid/Para.Beta, dataT[idx, :].real, "--",
                 label="$k/k_F={0:.4f}$, spectral".format(k/Para.kF))
        ax2.plot([TauGrid[0], TauGrid[-1]], [0.5, 0.5], linestyle='--')
        ax2.set_xlim([TauGrid[0]/Para.Beta, TauGrid[-1]/Para.Beta])
        ax2.set_xlabel("$\\tau/\\beta$")

    plt.legend(loc=1, frameon=False, fontsize=size)
    if Save:
        plt.savefig("DeltaG_T.pdf")
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

    # print abs(MomGrid-Para.kF)
    arr = np.amin(abs(MomGrid-Para.kF))
    kFidx = np.where(abs(arr - abs(MomGrid-Para.kF)) < 1.0e-20)[0][0]
    # print kFidx
    print "Mu=", Static[kFidx]
    Static -= Static[kFidx]  # subtract the self-energy shift

    print "Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr))

    print "MomGrid idx at the Fermi surface: ", kFidx
    # PlotSigmaW(Dynamic, MomGrid, kFidx, False)

    SigmaW, _ = Fourier.SpectralT2W(Dynamic)

    s0, s1 = SigmaW[kFidx, MaxFreq-1], SigmaW[kFidx, MaxFreq]
    print "Z=", 1.0-(s1.imag-s0.imag)/(2.0*np.pi/Para.Beta)
    dMu = (s0.real+s1.real)/2.0
    print "Dynamic chemical shift: ", dMu
    # PlotDataK(SigmaW, MaxFreq, Freq, False)

    BareG = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)

    for i, q in enumerate(MomGrid):
        BareG[i, :] = 1.0/(1j*phyFreq+q*q-Static[i]-Para.EF)

    # dG_W = BareG
    dG_W = 1.0/(1.0/BareG+SigmaW-dMu)
    # dG_W = (SigmaW-dMu)*BareG*BareG/(1-(SigmaW-dMu)*BareG)

    dG_T, _ = Fourier.SpectralW2T(dG_W)
    dG_Tp = Fourier.naiveW2T(dG_W)
    # print dG_T.shape
    # PlotSigmaT(dG_W, range(0, Para.MomGridSize, Para.MomGridSize/8), False)
    # print MomGrid[kFidx]/Para.kF
    # print BareG[kFidx, :]
    # PlotSigmaT(dG_W, [kFidx/2, kFidx, kFidx*2], False)
    PlotG(dG_T, kFidx, False)
    print "Maximum Error of \delta G: ", np.amax(abs(dG_T-dG_Tp))

    with open("dispersion.data", "w") as f:
        for k in range(Para.MomGridSize):
            f.write("{0} ".format(Static[k]))
        f.write("\n")

    with open("green.data", "w") as f:
        for k in range(Para.MomGridSize):
            for t in range(Para.TauGridSize):
                f.write("{0} ".format(dG_T[k, t].real))
        f.write("\n")

    ####### Tau space test  ############################
    # BareG = np.zeros((Para.MomGridSize, Para.TauGridSize), dtype=float)

    # for i, k in enumerate(MomGrid):
    #     Ek = k*k-Para.EF
    #     BareG[i, ] = np.exp(-Ek*TauGrid)/(1.0+np.exp(-Ek*Para.Beta))

    # PlotSigmaW(BareG, MomGrid, kFidx, False)
    # BareGW, _ = Fourier.SpectralT2W(BareG)
    # print BareGW.shape
    # BareGT, _ = Fourier.SpectralW2T(BareGW)
    # PlotSigmaT(BareGW, [kFidx-10, ], False)
    # PlotGT(BareGT, False)
