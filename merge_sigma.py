from scipy import integrate
from scipy import interpolate
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


def FourierT2W(SigmaT, FreqGrid, TauGrid, Beta):
    # print SigmaT.shape
    MomGridSize = SigmaT.shape[0]
    FreqGridSize = len(FreqGrid)
    # TauGridSize = len(TauGrid)

    SigmaT = np.array(SigmaT)

    C0 = SigmaT[:, 0]+SigmaT[:, -1]
    # C0 = SigmaT[:, 0]
    # C0 /= 1.2
    C1 = (SigmaT[:, 1]-SigmaT[:, 0])/(TauGrid[1]-TauGrid[0])
    C1 += (SigmaT[:, -1]-SigmaT[:, -2])/(TauGrid[-1]-TauGrid[-2])

    # for i in range(MomGridSize):
    #     SigmaT[i, :] -= C0[i]/2
    print SigmaT[12, 0], SigmaT[12, -1]
    # SigmaT[i, :] -= C1[i]*(Beta/2.0-TauGrid)/2.0

    SigmaW = np.zeros((MomGridSize, FreqGridSize), dtype=complex)

    for n, wn in enumerate(FreqGrid):
        tmp = SigmaT*np.exp(1j*TauGrid[np.newaxis, :]*wn)
        # SigmaW[:, n] = integrate.trapz(tmp, TauGrid, axis=1)
        SigmaW[:, n] = np.sum(tmp, axis=1)*Beta/len(TauGrid)
        # if wn > -29.5*np.pi/Beta and wn < 29.5*np.pi/Beta:
        #     SigmaW[:, n] = integrate.simps(tmp, TauGrid, axis=1)
        # else:
        #     SigmaW[:, n] = 0.0000/wn**2

        # SigmaW[:, n] = C1/wn**2

    return SigmaW, C0, C1


def FourierW2T(SigmaW, FreqGrid, TauGrid, Beta, C0, C1):
    # print SigmaT.shape
    MomGridSize = SigmaW.shape[0]
    TauGridSize = len(TauGrid)

    SigmaT = np.zeros((MomGridSize, TauGridSize), dtype=complex)

    for i, t in enumerate(TauGrid):
        Phase = np.exp(-1j*FreqGrid*t)
        SigmaT[:, i] = np.sum(SigmaW*Phase[np.newaxis, :], axis=1)/Beta

    # ExtFreqGrid = np.array(range(-10000, 20))
    # ExtFreqGrid = (2.0*ExtFreqGrid+1.0)*np.pi/Beta
    # for i, t in enumerate(TauGrid):
        # Phase = np.exp(-1j*FreqGrid*t)

    # for i in range(MomGridSize):
    #     SigmaT[i, :] += C0[i]/2
        # SigmaT[i, :] += C1[i]*(Beta/2.0-TauGrid)/2.0

    return SigmaT.real


if __name__ == "__main__":

    Para = param()
    Order = range(0, Para.Order+1)
    TauGrid = BuildTauGrid(Para, TauGridSize)
    MomGrid = BuildMomGrid(Para.MaxExtMom, MomGridSize)

    folder = "./Beta{0}_rs{1}_lambda{2}/".format(
        int(Para.Beta*Para.EF), Para.Rs, Para.Mass2)

    filename = "sigma_pid[0-9]+.dat"

    shape = (Para.Order+1, MomGridSize, TauGridSize)

    Data, Norm, Step = LoadFile(folder, filename, shape)

    Static, StaticErr = SigmaStatic(Data, Norm, Para)
    Dynamic, DynErr = SigmaT(Data, Norm, Para)

    # UniTauGrid = np.array(range(0, 1001), dtype=float)/1000.0*Para.Beta
    # UniTauGrid[0] += 2.0e-8
    # UniTauGrid[-1] -= 2.0e-8
    # f = interpolate.interp1d(TauGrid, Dynamic)
    # SigmaT = f(UniTauGrid)

    MaxFreq = 512
    Freq = np.array(range(-MaxFreq, MaxFreq))
    # print len(Freq)
    phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency

    # print SigmaT.shape
    SigmaW, C0, C1 = FourierT2W(Dynamic, phyFreq, TauGrid, Para.Beta)

    # print len(Freq)
    phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency
    SigmaT = FourierW2T(SigmaW, phyFreq, TauGrid, Para.Beta, C0, C1)

    # Plot SigmaW
    fig, ax = plt.subplots()
    idx = 12
    # ax.errorbar(TauGrid/Para.Beta, Dynamic[idx, :], yerr=0.0, fmt='o-',
    #             capthick=1, capsize=4, color=ColorList[idx], label="k={0}".format(MomGrid[idx]))

    # ax.errorbar(UniTauGrid/Para.Beta, SigmaT[idx, :], yerr=0.0, fmt='o-',
    #             capthick=1, capsize=4, color=ColorList[idx+1], label="k={0}, inter".format(MomGrid[idx]))

    # ax.errorbar(Freq, SigmaW[idx, :].imag, yerr=0.0, fmt='o-',
    #             capthick=1, capsize=4, color=ColorList[idx], label="k={0}".format(MomGrid[idx]))

    # ax.errorbar(Freq, -5.0/phyFreq**2, yerr=0.0, fmt='o-',
    #             capthick=1, capsize=4, color=ColorList[idx+1], label="k={0}".format(MomGrid[idx]))

    # ax.errorbar(Freq, C0[idx]/(phyFreq), yerr=0.0, fmt='o-',
    #             capthick=1, capsize=4, color=ColorList[idx+2], label="k={0}, 1/iwn".format(MomGrid[idx]))

    # ax.set_xlim([-1000, 1000])
    # ax.set_ylim([-0.1, 0.1])

    ax.errorbar(TauGrid/Para.Beta, Dynamic[idx, :], yerr=0.0, fmt='o-',
                capthick=1, capsize=4, color=ColorList[idx+1], label="k={0}, real".format(MomGrid[idx]))

    ax.errorbar(TauGrid/Para.Beta, SigmaT[idx, :], yerr=0.0, fmt='o-',
                capthick=1, capsize=4, color=ColorList[idx+2], label="k={0}, imag".format(MomGrid[idx]))

    # BareG = np.zeros((MomGridSize, 2*MaxFreq), dtype=complex)
    # for i, q in enumerate(MomGrid):
    #     BareG[i, :] = -1.0/(1j*phyFreq-q*q+Para.EF)

    # dG_W = SigmaW*BareG*BareG/(1-SigmaW*BareG)

    # dG_T = np.zeros((MomGridSize, TauGridSize), dtype=float)

    # ax.set_xticks([0.0,0.04,0.08,0.12])
    # ax.set_yticks([0.35,0.4,0.45,0.5])
    # ax.set_ylim([-0.02, 0.125])
    # ax.set_ylim([0.07, 0.125])
    # ax.xaxis.set_label_coords(0.97, -0.01)
    # # ax.yaxis.set_label_coords(0.97, -0.01)
    # ax.text(-0.012,0.52, "$-I$", fontsize=size)
    # ax.set_ylabel("$-\Gamma_4(\omega=0, q)$", size=size)

    # ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
# plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
