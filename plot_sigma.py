#!/usr/bin/env python3
from utility.IO import *
import utility.fourier as fourier
from utility.plot import *
import argparse

parser = argparse.ArgumentParser("Specify some parameters.")
parser.add_argument("folder")
args = parser.parse_args()

folder = args.folder
print("Folder to plot : " + folder)


def PlotSigmaW_RI(SigmaT, MomGrid, idx, Save=True):
    SigmaW, Spec = Fourier.SpectralT2W(SigmaT[idx, :])

    _, (ax1, ax2) = plt.subplots(1, 2)
    Nw = int(len(phyFreq)/4)


    ax1.plot(phyFreq[Nw:-Nw], SigmaW.real[Nw:-Nw], "r-", label="real")
    
    # ax1.set_ylim([-0.05, 0.05])
    # ax2.set_ylim([1.0e-6, 1.1])
    # ax2.set_yscale("log")
    ax2.plot(phyFreq[Nw:-Nw], SigmaW.imag[Nw:-Nw], "b-", label="imaginary")

    ax1.set_ylabel("$\Sigma(i\\omega_n, k_F)$")
    ax1.set_xlabel("$\\omega_n$")
    # ax2.set_ylabel("$\\Sigma(\\tau, k)$")
    ax2.set_ylabel("$\Sigma(i\\omega_n, k_F)$")
    ax2.set_xlabel("$\\omega_n$")
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax1.legend(loc=1, frameon=False, fontsize=size)
    ax2.legend(loc=1, frameon=False, fontsize=size)
    plt.tight_layout()
    plt.show()


def PlotSigmaT_RI(SigmaT, MomGrid, idx, Save=True):
    _, (ax1, ax2) = plt.subplots(1, 2)
    Nw = int(len(phyFreq)/4)

    ax1.plot(TauGrid/Para.Beta, SigmaT.real[idx, :], "r-", label="real")
    
    # ax1.set_ylim([-0.05, 0.05])
    # ax2.set_ylim([1.0e-6, 1.1])
    # ax2.set_yscale("log")
    ax2.plot(TauGrid/Para.Beta, SigmaT.imag[idx, :], "b-", label="imaginary")

    ax1.set_ylabel("$\Sigma(\\tau, k_F)$")
    ax1.set_xlabel("$\\tau/\\beta$")
    # ax2.set_ylabel("$\\Sigma(\\tau, k)$")
    ax2.set_ylabel("$\Sigma(\\tau, k_F)$")
    ax2.set_xlabel("$\\tau/\\beta$")
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax1.legend(loc=1, frameon=False, fontsize=size)
    ax2.legend(loc=1, frameon=False, fontsize=size)
    plt.tight_layout()
    plt.show()



def PlotSigmaW(SigmaT, MomGrid, idx, Save=True):
    SigmaW, Spec = Fourier.SpectralT2W(SigmaT[idx, :])
    C0 = SigmaT[idx, 0]+SigmaT[idx, -1]
    SigmaWp = Fourier.naiveT2W(SigmaT[idx, :])

    _, (ax1, ax2) = plt.subplots(1, 2)
    # ax.plot(Freq, Spec[idx, :], "k--",
    #         label="k={0}, spectral".format(Grid.MomGrid[idx]))
    k = MomGrid[idx]/Para.kF
    print(phyFreq)
    ax1.plot(phyFreq, SigmaW.real, "r-",
             label="Spectral Fourier")
    ax1.plot(phyFreq, SigmaWp.imag, "b--",
             label="Naive Fourier")
    Nw = len(phyFreq)
    ax1.plot(phyFreq[:int(Nw*0.45)], -C0/phyFreq[:int(Nw*0.45)], "g--",
             label="Tail $\\sim 1/\omega_n$")
    ax1.plot(phyFreq[int(Nw*0.55):], -C0/phyFreq[int(Nw*0.55):], "g--")
    # ax1.plot(phyFreq, (-1/(1.0j*phyFreq)).imag, "g--",
    #          label="$k/k_F={0:.4f}$, tail".format(k))

    # ax1.plot(phyFreq, SigmaW.real, "bo",
    #          label="$k/k_F={0}$, spectral".format(k))
    # ax1.plot(phyFreq, SigmaWp.real, "b--",
    #          label="$k/k_F={0}$, naive".format(k))

    ax1.set_ylim([-0.05, 0.05])
    # ax1.set_ylim([-0.02, 0.02])
    ax2.set_ylim([1.0e-6, 1.1])
    # ax2.set_ylim([0.495, 0.505])
    ax2.set_yscale("log")

    SigmaTpn = Fourier.naiveW2T(SigmaW)
    SigmaTp, Spec = Fourier.SpectralW2T(SigmaW)
    ax2.plot(TauGrid/Para.Beta, SigmaT[idx, :], "ko", markersize=2,
             label="Exact")
    ax2.plot(TauGrid/Para.Beta, SigmaTpn.real, "b--", lw=1,
             label="Naive Fourier")
    ax2.plot(TauGrid/Para.Beta, SigmaTp.real, "r-",
             label="Spectral Fourier")
    # ax2.set_xlim([0.0-1.0e-3, 1.0+1.0e-3])

    # ax1.set_ylabel("$\\Sigma(i\\omega_n, k)$")
    ax1.set_ylabel("$G(i\\omega_n, k)$")
    ax1.set_xlabel("$\\omega_n$")
    # ax2.set_ylabel("$\\Sigma(\\tau, k)$")
    ax2.set_ylabel("$G(\\tau, k)$")
    ax2.set_xlabel("$\\tau/\\beta$")
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax1.legend(loc=1, frameon=False, fontsize=size)
    ax2.legend(loc=1, frameon=False, fontsize=size)
    plt.tight_layout()
    plt.savefig(f"sigma_k={k}.pdf")
    # plt.grid()
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
    plt.grid()
    _ = InteractiveLegend(ax)
    if Save:
        plt.savefig("DataK.pdf")
    else:
        plt.show()


def PlotSigmaT(dataMomTau, kList, Save=True):
    fig, ax = plt.subplots()
    dataW, _ = Fourier.SpectralT2W(dataMomTau)
    # dataWp = Fourier.naiveT2W(dataT)

    for idx in kList:
        idx = int(idx)
        k = MomGrid[idx]
        dataT = dataMomTau[idx]
        # dataT, _ = Fourier.SpectralW2T(dataW[idx, :])
        dataTp = Fourier.naiveW2T(dataW[idx, :])
        ax.plot(TauGrid/Para.Beta, dataT.real, "r-",
                label="original: $k={:.2f}k_F$".format(k/Para.kF))
        ax.plot(TauGrid/Para.Beta, dataTp.real, "b--",
                label="Naive Fourier: $k={:.2f}k_F$".format(k/Para.kF))
        Ek = k*k-Para.EF
        
        # G0 = np.exp(-Ek*TauGrid)/(1.0+np.exp(-Ek*Para.Beta))
        # ax.plot(TauGrid/Para.Beta, G0,
                # "ko", markersize=1, label="Bare G: k={:.2f}".format(k/Para.kF))

    ax.set_xlim([TauGrid[0]/Para.Beta, TauGrid[-1]/Para.Beta])
    ax.set_xlabel("$\\tau/\\beta$")
    ax.set_ylabel("$\\Sigma(\\tau)$")

    plt.legend(loc=1, frameon=False, fontsize=size)
    # plt.grid()
    _ = InteractiveLegend(ax)
    if Save:
        plt.savefig("DeltaG.pdf")
    else:
        plt.show()


def PlotG(dataT, kFidx, Save=True):
    _, (ax1, ax2) = plt.subplots(1, 2)


    ax1.plot(MomGrid/Para.kF, dataT[:, -1], "ro-",
             label="$\\tau=0^-$")
    ax1.plot(MomGrid/Para.kF, dataT[:, 0], "bs-",
             label="$\\tau=0^+$")
    # ax.plot(TauGrid, dataTp.real, "--",
    #         label="$k/k_F={0}$, naive".format(k))

    ax1.plot([MomGrid[0], MomGrid[-1]], [0.5, 0.5], linestyle='--')
    ax1.axvline(x=1.0, linestyle='--')
    ax1.set_xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    ax1.set_xlabel("$k/k_F$")
    # ax1.set_ylim([-0.12, 0.12])

    for idx in [kFidx/2, kFidx, kFidx+kFidx/2]:
        k = MomGrid[int(idx)]
        ax2.plot(TauGrid/Para.Beta, dataT[int(idx), :].real, "--",
                 label="$k/k_F={0:.2f}$, spectral".format(k/Para.kF))
        ax2.set_xlim([TauGrid[0]/Para.Beta, TauGrid[-1]/Para.Beta])
        ax2.set_xlabel("$\\tau/\\beta$")
    ax2.plot([TauGrid[0], TauGrid[-1]], [0.5, 0.5], linestyle='--')

    plt.legend(loc=1, frameon=False, fontsize=size)
    plt.grid()
    if Save:
        plt.savefig("DeltaG_T.pdf")
    else:
        plt.show()


def PlotGT(dataW, Save=True):
    _, ax = plt.subplots()

    dataT, _ = Fourier.SpectralW2T(dataW)
    dataTp = Fourier.naiveW2T(dataW)

    ax.plot(MomGrid/Para.kF, dataT[:, -1].real, "r-",
            label="Spectral Fourier, $\\tau=-0^+$")
    ax.plot(MomGrid/Para.kF, dataTp[:, -1].real, "b--",
            label="Naive Fourier, $\\tau=-0^+$")

    Ek = MomGrid**2-Para.EF
    # print "Ek: ", Ek
    # G0 = np.exp(-Ek*TauGrid)/(1.0+np.exp(-Ek*Para.Beta))
    n0 = 1.0/(1.0+np.exp(Ek*Para.Beta))
    ax.plot(MomGrid/Para.kF, n0,
            "ko", markersize=2, label="Exact")
    # ax1.plot(MomGrid/Para.kF, dataT[:, 0], "bs-",
    #          label="$\\tau=0^+$")
    # ax.plot(TauGrid, dataTp.real, "--",
    #         label="$k/k_F={0}$, naive".format(k))

    # ax.plot([MomGrid[0], MomGrid[-1]], [0.5, 0.5], linestyle='--')
    # ax.axvline(x=1.0, linestyle='--')
    ax.set_xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    ax.set_ylabel("$n(k)$")
    ax.set_xlabel("$k/k_F$")
    # ax1.set_ylim([-0.12, 0.12])

    plt.legend(loc=1, frameon=False, fontsize=size)
    # plt.grid()
    plt.tight_layout()
    plt.savefig("FermiDist.pdf")
    if Save:
        plt.savefig("DeltaG_T.pdf")
    else:
        plt.show()


def check_convergence(Para, Data):
    print(Para.Order)
    Order = range(0, Para.Order+1)
    print("    Q/kF,      Data,      Error")
    for o in Order[1:]:
        print("Order {0}".format(o))
        # sum all orders
        DataAllList = [np.sum(d[1:o+1, ...], axis=0) for d in Data]
        # average tau
        DataAllList = [np.average(d[:,:], axis=1) for d in DataAllList]
        
        Dat, Err = Estimate(DataAllList, Norm)
        print("tau average Sigma({0:2.1f}kF):, {1:10.6f}, {2:10.6f}".format(
            MomGrid[0], Dat[0], Err[0]))
    # sys.exit(0)



if __name__ == "__main__":
    while True:
        Para = param(folder)
        Order = range(0, Para.Order+1)

        MaxFreq = 3000
        Freq = np.array(range(-MaxFreq, MaxFreq))
        phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency

        shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)
        Data, Norm, Step, Grids = LoadFile(folder, "sigma_pid[0-9]+.dat", shape)

        # Data.shape : (pid, order+1, MomGrid, TauGrid)
        
        TauGrid = Grids["TauGrid"]
        MomGrid = Grids["KGrid"]

        check_convergence(Para, Data)

        Fourier = fourier.fourier(TauGrid, phyFreq, Para.Beta)
        Fourier.InitializeKernel(100.0, 1024, "Fermi", 1.0e-13)

        # first order is a constant of tau
        # Static(MomGrid) : order 1
        Static, StaticErr = Estimate(
            Data, Norm, lambda d: np.average(d[1, :, :], axis=1))
        # Dynamic(MomGrid, TauGrid) : order 2-n
        Dynamic, DynErr = Estimate(
            Data, Norm, lambda d: np.sum(d[1:Para.Order+1, ...], axis=0))
        # Static_MomTau(MomGrid, TauGrid) : order 2-n
        # Static_MomTau, StaticErr = Estimate(Data, Norm, lambda d: d[1, :, :])

        
        arr = np.amin(abs(MomGrid-Para.kF))
        kFidx = np.where(abs(arr - abs(MomGrid-Para.kF)) < 1.0e-20)[0][0]
        
        print("Mu=", Static[kFidx])
        Static -= Static[kFidx]  # subtract the self-energy shift

        print("Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr)))

        print("MomGrid idx at the Fermi surface:{0}, KF: {1}=={2}".format(kFidx,MomGrid[kFidx],Para.kF))
        # PlotSigmaW(Dynamic, MomGrid, kFidx, False)

        SigmaW, _ = Fourier.SpectralT2W(Dynamic)

        s0, s1 = SigmaW[kFidx, MaxFreq-1], SigmaW[kFidx, MaxFreq]
        print("Z=", 1.0-(s1.imag-s0.imag)/(2.0*np.pi/Para.Beta))
        dMu = (s0.real+s1.real)/2.0
        print("Dynamic chemical shift: ", dMu)

        BareG = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)

        for i, q in enumerate(MomGrid):
            BareG[i, :] = 1.0/(1j*phyFreq+q*q-Static[i]-Para.EF)

        # dG_W = BareG
        dG_W = 1.0/(1.0/BareG+SigmaW-dMu)
        # dG_W = (SigmaW-dMu)*BareG*BareG/(1-(SigmaW-dMu)*BareG)

        dG_T, _ = Fourier.SpectralW2T(dG_W)
        dG_Tp = Fourier.naiveW2T(dG_W)
        # PlotSigmaT(dG_W, range(0, Para.MomGridSize, Para.MomGridSize/8), False)
        # print MomGrid[kFidx]/Para.kF
        # print BareG[kFidx, :]
        # PlotSigmaT(Static_MomTau, [kFidx], False)
        # PlotG(dG_T, kFidx, False)
        print("Maximum Error of \delta G: ", np.amax(abs(dG_T-dG_Tp)))

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
        BareG = np.zeros((Para.MomGridSize, Para.TauGridSize), dtype=float)

        for i, k in enumerate(MomGrid):
            Ek = k*k-Para.EF
            BareG[i, ] = np.exp(-Ek*TauGrid)/(1.0+np.exp(-Ek*Para.Beta))

        PlotSigmaT_RI(Dynamic, MomGrid, kFidx, False)
        # BareGW, _ = Fourier.SpectralT2W(BareG)
        # BareGT, _ = Fourier.SpectralW2T(BareGW)
        
        break
