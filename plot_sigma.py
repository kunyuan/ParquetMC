#!/usr/bin/env python3
from utility.IO import *
import utility.fourier as fourier
from utility.plot import *
import argparse
from mpl_toolkits.mplot3d import Axes3D


parser = argparse.ArgumentParser("Specify some parameters.")
parser.add_argument("folder")
args = parser.parse_args()

folder = args.folder
print("Folder to plot : " + folder)

def AnaliticFock(k, l):
    kF = Para.kF
    # Analytic Fock energy
    y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((k-kF)/l)-l/kF*np.arctan((k+kF)/l) -
                      (l*l-k*k+kF*kF)/4.0/k/kF*np.log((l*l+(k-kF)**2)/(l*l+(k+kF)**2)))

    k = kF
    Mu = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((k-kF)/l)-l/kF*np.arctan((k+kF)/l) -
                       (l*l-k*k+kF*kF)/4.0/k/kF*np.log((l*l+(k-kF)**2)/(l*l+(k+kF)**2)))
    return y-Mu


def PlotStatic():
    fig, ax1 = plt.subplots()
    fock = [AnaliticFock(k, Para.Lambda) for k in MomGrid]

    ax1.plot(MomGrid/Para.kF, fock, "r-", label="Analitical Fock")
    ax1.plot(MomGrid/Para.kF, Static, "b-", label="Static")

    plt.legend(loc=1, frameon=False, fontsize=size)
    plt.show()


def CoulombFock(k, l):
    kF = Para.kF
    y = 2.0*kF/np.pi*(1.0 + (k*k-kF*kF)/(2*k*kF)*np.log(np.abs((k-kF)/(k+kF))))
    return y

def PlotSigmaW_RI(SigmaT, MomGrid, idx, Save=False):
    SigmaW, Spec = Fourier.SpectralT2W(SigmaT[idx, :])

    # Fourier1 = fourier.fourier(dataTau[:,0], phyFreq, Para.Beta)
    # Fourier1.InitializeKernel(100.0, 1024, "Fermi", 1.0e-13)

    # dataWFromTau, _ = Fourier1.SpectralT2W(dataTau[:,1])

    _, (ax1, ax2) = plt.subplots(1, 2)
    Nw = int(len(phyFreq)/2) - 94


    ax1.plot(phyFreq[Nw:-Nw]/Para.EF, SigmaW.real[Nw:-Nw], "r-", label="real")
    
    ax2.plot(phyFreq[Nw:-Nw]/Para.EF, SigmaW.imag[Nw:-Nw], "r-", label="imaginary")

    print(np.max(SigmaW.real),np.max(SigmaW.imag),np.max(SigmaW.real)/np.max(SigmaW.imag))
    print((np.max(SigmaW.real)-np.min(SigmaW.real)) , (np.max(SigmaW.imag)-np.min(SigmaW.imag)),
     (np.max(SigmaW.real)-np.min(SigmaW.real)) / (np.max(SigmaW.imag)-np.min(SigmaW.imag)) )

    ax1.set_ylabel("$\Sigma(i\\omega_n, k_F)$")
    ax1.set_xlabel("$\\omega_n/E_F$")
    # ax2.set_ylabel("$\\Sigma(\\tau, k)$")
    ax2.set_ylabel("$\Sigma(i\\omega_n, k_F)$")
    ax2.set_xlabel("$\\omega_n/E_F$")
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax1.legend(loc=1, frameon=False, fontsize=size)
    ax2.legend(loc=1, frameon=False, fontsize=size)
    plt.tight_layout()
    plt.show()


def PlotSigmaT_RI(SigmaT, MomGrid, idx, Save=True):
    _, (ax1, ax2) = plt.subplots(2)
    Nw = int(len(phyFreq)/4)

    print(SigmaT.real[idx, :])
    # print(dataTau[:,0]/Para.Beta, 0.50087*dataTau[:,1])

    # with open("taugrid.data", "w") as ff:
    #     for i in range(0, len(TauGrid)):
    #         ff.write(str(TauGrid[i]) + "    " +str(SigmaT.real[idx, i]) + "\n")

    ax1.plot(TauGrid/Para.Beta, SigmaT.real[idx, :], "r-", label="real")
    ax2.plot(TauGrid/Para.Beta, SigmaT.imag[idx, :], "r-", label="real")
    # ax1.plot(dataTau[:,0]/Para.Beta, 0.5*dataTau[:,1], "b-", label="Wang Tao")
    
    ax1.set_ylabel("$\Sigma(\\tau, k_F)$")
    ax1.set_xlabel("$\\tau/\\beta$")
    ax1.legend(loc=1, frameon=False, fontsize=size)

    ax2.set_ylabel("$\Sigma(\\tau, k_F)$")
    ax2.set_xlabel("$\\tau/\\beta$")
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
             label="Ekact")
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
    # plt.savefig(f"sigma_k={k}.pdf")
    # plt.grid()
    if Save:
        plt.savefig("SigmaW.pdf")
    else:
        plt.show()


def PlotDatax(datax, wn, FreqGrid, Save=True):
    fig, ax = plt.subplots()
    ax.plot(MomGrid, datax[:, wn].real, "ro-",
            label="$n={0}$, real".format(FreqGrid[wn]))
    ax.plot(MomGrid, datax[:, wn].imag, "bo-",
            label="$n={0}$, imag".format(FreqGrid[wn]))

    plt.axvline(k=Para.kF, linestyle='--')

    plt.legend(loc=1, frameon=False, fontsize=size)
    plt.grid()
    _ = InteractiveLegend(ax)
    if Save:
        plt.savefig("Datax.pdf")
    else:
        plt.show()


def PlotSigmaT(dataMomTau, kList, Save=False):
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
        # ax.plot(TauGrid/Para.Beta, dataTp.real, "b--",
                # label="Naive Fourier: $k={:.2f}k_F$".format(k/Para.kF))
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

def CheckFock():
    fig, ax1 = plt.subplots()
    MomGrid = np.array([1.0*i/100.0*Para.kF for i in range(1, 1700)])
    # fock = [AnaliticFock(k, np.sqrt(Para.Mass2 + Para.Lambda)) for k in MomGrid]
    fock = [AnaliticFock(k, 0.00000000001) for k in MomGrid]
    fockCoulomb = [CoulombFock(k, 0.0000000001) for k in MomGrid]

    # ax1.plot(MomGrid/Para.kF, Static, "k-", label="Fock in code")
    ax1.plot(MomGrid/Para.kF, fock, "c-", label="Analitical Fock")
    # ax1.plot(MomGrid/Para.kF, fockCoulomb, "r-", label="Coulomb Fock")
    plt.legend(loc=1, frameon=False, fontsize=size)

    plt.show()


def GTPlot(kList, Save=False):
    fig, ax1 = plt.subplots()
    BareGTFW, _ = Fourier.SpectralW2T(BareGW)
    BareGT = np.zeros((Para.MomGridSize, len(TauGrid)), dtype=complex)
    for i, k in enumerate(MomGrid):
            Ek = k*k-Para.EF
            BareGT[i, ] = np.exp(-Ek*TauGrid)/(1.0+np.exp(-Ek*Para.Beta))
    # fig = plt.figure()
    # ax = axes3D(fig)
    # k, Y = np.meshgrid(TauGrid, MomGrid)

    # ax.plot_surface(k, Y, dG_T.real, rstride=1, cstride=1, cmap='rainbow')
    
    # with open(os.path.join(folder, "green2.data"), "w") as f:
    #         for k in range(Para.MomGridSize):
    #             for t in range(Para.TauGridSize):
    #                 f.write("{0} ".format(dG_T[k, t].real))
    #         f.write("\n")

    # d1 = np.loadtxt(os.path.join(folder, "green.data"), dtype=float).reshape(64,128)

    # d1 = np.loadtxt(os.path.join(folder,"selfconsistent/green_order2.data"), dtype=float).reshape(64,128)
    d2 = np.loadtxt(os.path.join(folder,"selfconsistent/green_order3.data"), dtype=float).reshape(64,128)
    # d3 = d1 + d2

    for idx in kList:
        idx = int(idx)
        k = MomGrid[idx]
        ax1.plot(TauGrid/Para.Beta, dG_T[idx].real, "k-", label="dG")
        # ax1.plot(TauGrid/Para.Beta, d1[idx].real, "r-", label="d1")
        ax1.plot(TauGrid/Para.Beta, d2[idx].real, "r-", label="d2")
        # ax1.plot(TauGrid/Para.Beta, d3[idx].real, "c-", label="d3")

        # ax1.plot(TauGrid/Para.Beta, BareGTFW[idx], "b-", label="BareG From W")
        # ax1.plot(TauGrid/Para.Beta, BoldG_T[idx], "r-", label="BoldG")
        # ax1.plot(TauGrid/Para.Beta, BareGT[idx], "c-", label="GBare")
        
        # ax1.plot(TauGrid/Para.Beta, BareG_T[idx], "b-", label="GBare")
    plt.legend(loc=1, frameon=False, fontsize=size)

    if Save:
        plt.savefig("GTPlot.pdf")
    else:
        plt.show()


def PlotG(dataMomTau, kList, Save=False):
    _, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(MomGrid/Para.kF, dataMomTau[:, -1], "r-",
             label="$\\tau=0^-$")
    ax1.plot(MomGrid/Para.kF, dataMomTau[:, 0], "b-",
             label="$\\tau=0^+$")
    # ax.plot(TauGrid, dataTp.real, "--",
    #         label="$k/k_F={0}$, naive".format(k))

    ax1.plot([MomGrid[0], MomGrid[-1]], [0.5, 0.5], linestyle='--')
    ax1.axvline(k=1.0, linestyle='--')
    ax1.set_xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    ax1.set_xlabel("$k/k_F$")
    ax1.set_ylabel("$G(k,\\tau=0^{\\pm})$")
    # ax1.set_ylim([-0.12, 0.12])

    for idx in kList:
        k = MomGrid[int(idx)]
        ax2.plot(TauGrid/Para.Beta, dataMomTau[int(idx), :].real, "--",
                 label="$k/k_F={0:.2f}$, spectral".format(k/Para.kF))
        ax2.set_xlim([TauGrid[0]/Para.Beta, TauGrid[-1]/Para.Beta])
    ax2.set_xlabel("$\\tau/\\beta$")
    ax2.set_ylabel("$G(\\tau)$")
    ax2.plot([TauGrid[0], TauGrid[-1]], [0.5, 0.5], linestyle='--')

    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()

    plt.legend(loc=1, frameon=False, fontsize=size)
    plt.grid()
    if Save:
        plt.savefig("DeltaG_T.pdf")
    else:
        plt.show()


def PlotGT(dataW, Save=False):
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
            "ko", markersize=2, label="Ekact")
    # ax1.plot(MomGrid/Para.kF, dataT[:, 0], "bs-",
    #          label="$\\tau=0^+$")
    # ax.plot(TauGrid, dataTp.real, "--",
    #         label="$k/k_F={0}$, naive".format(k))

    # ax.plot([MomGrid[0], MomGrid[-1]], [0.5, 0.5], linestyle='--')
    # ax.axvline(k=1.0, linestyle='--')
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
    # sys.ekit(0)



if __name__ == "__main__":
    while True:

        # dataFreqReal = np.loadtxt(os.path.join(folder, "sigma_freq_real.dat"))
        # dataFreqImag = np.loadtxt(os.path.join(folder, "sigma_freq_imag.dat"))
        # dataTau = np.loadtxt(os.path.join(folder, "sigma_tau.dat"))

        # with open("taugrid.data", "w") as ff:
        #     for i in range(0, len(dataTau[:,0])):
        #         ff.write(str(dataTau[i,0]) + "\n")

        # sys.ekit(0)

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

        # fig, ax = plt.subplots()
        # ax.plot(range(len(TauGrid)), TauGrid, 'r-')
        # ax.plot(range(len(dataTau[:,0])), dataTau[:,0], 'b-', label="Wang Tao")
        # ax.legend(loc=1, frameon=False, fontsize=size)
        # plt.tight_layout()
        # plt.show()

        # # check_convergence(Para, Data)
        # sys.ekit()

        

        Fourier = fourier.fourier(TauGrid, phyFreq, Para.Beta)
        Fourier.InitializeKernel(100.0, 1024, "Fermi", 1.0e-13)

        # first order is a constant of tau
        # Static(MomGrid) : order 1
        Static, StaticErr = Estimate(
            Data, Norm, lambda d: np.average(d[1, :, :], axis=1))

        # Dynamic(MomGrid, TauGrid) : order 2-n
        Dynamic, DynErr = Estimate(
            Data, Norm, lambda d: np.sum(d[2:3, ...], axis=0))
        # Dynamic, DynErr = Estimate(
            # Data, Norm, lambda d: d[2, ...])
        # Static_MomTau(MomGrid, TauGrid) : order 2-n
        Static_MomTau, StaticErr = Estimate(Data, Norm, lambda d: d[1, :, :])
        
        arr = np.amin(abs(MomGrid-Para.kF))
        kFidx = np.where(abs(arr - abs(MomGrid-Para.kF)) < 1.0e-20)[0][0]
        
        SigmaW, _ = Fourier.SpectralT2W(Static_MomTau)
        s0, s1 = SigmaW[kFidx, MaxFreq-1], SigmaW[kFidx, MaxFreq]
        Z = 1.0-(s1.imag-s0.imag)/(2.0*np.pi/Para.Beta)
        print("Z=", Z)

        # PlotSigmaT_RI(-Dynamic, MomGrid, kFidx, Save=False)
        # PlotSigmaW_RI(-Static_MomTau, MomGrid, kFidx, Save=False)
        # sys.ekit(0)


        print("Mu=", Static[kFidx])
        Static -= Static[kFidx]  # subtract the self-energy shift

        PlotStatic()
        sys.exit()

        print("Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr)))

        print("MomGrid idx at the Fermi surface:{0}, KF: {1}=={2}".format(kFidx,MomGrid[kFidx],Para.kF))
        # PlotSigmaW(Dynamic, MomGrid, kFidx, False)

        allorder, DynErr = Estimate(
            Data, Norm, lambda d: np.sum(d[1:Para.Order+1, ...], axis=0))
        SigmaW, _ = Fourier.SpectralT2W(Dynamic)

        

        s0, s1 = SigmaW[kFidx, MaxFreq-1], SigmaW[kFidx, MaxFreq]
        Z = 1.0-(s1.imag-s0.imag)/(2.0*np.pi/Para.Beta)
        print("Z=", Z)
        dMu = (s0.real+s1.real)/2.0
        print("Dynamic chemical shift: ", dMu)

        BareGW = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)

        for i, q in enumerate(MomGrid):
            BareGW[i, :] = Z/(1j*phyFreq + (q*q-Para.EF) + Static[i])


        BoldGW = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)
        for i, q in enumerate(MomGrid):
            for j, w in enumerate(phyFreq):
                BoldGW[i, j] = Z/( 1j*w + (q*q-Para.EF) + Static[i] + SigmaW[i,j] )

        

        # CheckFock()
        # sys.ekit(0)

        dG_W = BoldGW - BareGW

        dG_T, _ = Fourier.SpectralW2T(dG_W)

        

        BoldG_T, _ = Fourier.SpectralW2T(BoldGW)
        BareG_T, _ = Fourier.SpectralW2T(BareGW)

        print(kFidx/2)
        # with open(os.path.join(folder, "green.data"), "w") as f:
        #     for k in range(Para.MomGridSize):
        #         for t in range(Para.TauGridSize):
        #             f.write("{0} ".format(dG_T[k, t].real))
        #     f.write("\n")

        # GTPlot([kFidx/2], False)
        

        PlotSigmaT(Dynamic, [kFidx/2, kFidx, kFidx+kFidx/2], False)


        print("dG_T shape:{0}".format(np.array(dG_T).shape))
        
        print("Maximum Error of \delta G: ", np.amax(abs(dG_T-dG_Tp)))

        # with open("dispersion.data", "w") as f:
        #     for k in range(Para.MomGridSize):
        #         f.write("{0} ".format(Static[k]))
        #     f.write("\n")

        

        ####### Tau space test  ############################
        BareG = np.zeros((Para.MomGridSize, Para.TauGridSize), dtype=float)

        for i, k in enumerate(MomGrid):
            Ek = k*k-Para.EF
            BareG[i, ] = np.exp(-Ek*TauGrid)/(1.0+np.exp(-Ek*Para.Beta))

        # PlotSigmaW(allorder, MomGrid, kFidx, False)
        # BareGW, _ = Fourier.SpectralT2W(BareG)
        # BareGT, _ = Fourier.SpectralW2T(BareGW)
        
        break
