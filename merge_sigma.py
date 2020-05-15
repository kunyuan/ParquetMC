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

    Coef = np.dot(Dynamic, Tau2Coef)
    FittedDyn = np.dot(Coef, Coef2Tau)

    print "Max of |Sigma-Fitted Sigma|: ", np.amax(abs(Dynamic-FittedDyn))

    # fig, ax = plt.subplots()
    # idx = 30
    # ax.errorbar(Grid.TauGrid/Para.Beta, Dynamic[idx, :], yerr=DynErr[idx, :], fmt='o-',
    #             capthick=1, capsize=4, color='r', label="k={0}, real".format(Grid.MomGrid[idx]))

    # ax.plot(Grid.TauGrid/Para.Beta,
    #         FittedDyn[idx, :], 'bo-', label="k={0}, real".format(Grid.MomGrid[idx]))

    # plt.legend(loc=1, frameon=False, fontsize=size)
    # plt.show()

    return Coef


if __name__ == "__main__":

    fig, ax = plt.subplots()

    Para = param()
    Grid = grid(Para)
    Spectral = spectral(Para)
    Order = range(0, Para.Order+1)

    folder = "./Data"

    filename = "sigma_pid[0-9]+.dat"

    shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)

    Data, Norm, Step = LoadFile(folder, filename, shape)

    Static, StaticErr = SigmaStatic(Data, Norm, Para)
    Dynamic, DynErr = SigmaT(Data, Norm, Para)

    print "Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr))

    ut, st, vt = Spectral.TauKernel(Grid.TauGrid, "Fermi")
    FittedDyn, Coef = FitData(Dynamic, Para.TauBasisSize, vt)
    spectral = Data2Spectral(Dynamic, Para.TauBasisSize, ut, st, vt)

    # idx = 30
    # ax.errorbar(Grid.TauGrid/Para.Beta, Dynamic[idx, :], yerr=DynErr[idx, :], fmt='o-',
    #             capthick=1, capsize=4, color='r', label="k={0}, real".format(Grid.MomGrid[idx]))

    # ax.plot(Grid.TauGrid/Para.Beta,
    #         FittedDyn[idx, :], 'bo-', label="k={0}, real".format(Grid.MomGrid[idx]))

    # plt.legend(loc=1, frameon=False, fontsize=size)
    # plt.show()

    MaxFreq = 512
    Freq = np.array(range(-MaxFreq, MaxFreq))
    # # print len(Freq)
    phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency

    uw, sw, vw = Spectral.MatFreqKernel(phyFreq, "Fermi")
    SigmaW = Spectral2Data(spectral, Para.TauBasisSize, uw, sw, vw)

    # idx = 12

    # ax.plot(Freq, SigmaW[idx, :].real, "ro-",
    #         label="k={0}, real".format(Grid.MomGrid[idx]))
    # ax.plot(Freq, SigmaW[idx, :].imag, "bs-",
    #         label="k={0}, imag".format(Grid.MomGrid[idx]))

    # plt.legend(loc=1, frameon=False, fontsize=size)
    # plt.show()

    BareG = np.zeros((Grid.MomSize, len(phyFreq)), dtype=complex)
    for i, q in enumerate(Grid.MomGrid):
        BareG[i, :] = -1.0/(1j*phyFreq-q*q+Para.EF)

    # dG_W = 1.0/(1.0/BareG-SigmaW)
    dG_W = SigmaW*BareG*BareG/(1-SigmaW*BareG)

    FittedGW, Coefp = FitData(dG_W, Para.TauBasisSize, vw)

    # idx = 30
    # ax.plot(Freq,
    #         dG_W[idx, :], 'rs-', label="k={0}, real".format(Grid.MomGrid[idx]))
    # ax.plot(Freq,
    #         FittedGW[idx, :], 'bo-', label="k={0}, real".format(Grid.MomGrid[idx]))

    # plt.legend(loc=1, frameon=False, fontsize=size)
    # plt.show()

    spectral = Data2Spectral(dG_W, Para.TauBasisSize, uw, sw, vw)

    print "Max imaginary part of the spectral: ", np.amax(abs(spectral.imag))
    dG_T = Spectral2Data(spectral, Para.TauBasisSize, ut, st, vt)

    idx = 32
    ax.plot(Grid.TauGrid/Para.Beta, dG_T[idx, :].real, "ro-",
            label="k={0}, real".format(Grid.MomGrid[idx]))

    plt.legend(loc=1, frameon=False, fontsize=size)
    plt.show()

    # ax.set_xticks([0.0,0.04,0.08,0.12])
    # ax.set_yticks([0.35,0.4,0.45,0.5])
    # ax.set_ylim([-0.02, 0.125])
    # ax.set_ylim([0.07, 0.125])
    # ax.xaxis.set_label_coords(0.97, -0.01)
    # # ax.yaxis.set_label_coords(0.97, -0.01)
    # ax.text(-0.012,0.52, "$-I$", fontsize=size)
    # ax.set_ylabel("$-\Gamma_4(\omega=0, q)$", size=size)

    # ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)
