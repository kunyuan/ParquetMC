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


# XType = "Tau"
XType = "Mom"
# XType = "Z"
# XType = "Freq"
OrderByOrder = False
# 0: I, 1: T, 2: U, 3: S

Para = param(folder)
Order = range(2, Para.Order+1)

Data, Norm, Step, Grid = LoadFile(folder, "sigma_pid[0-9]+.dat")

MomGrid = Grid["KGrid"]
TauGrid = Grid["TauGrid"]
MomGridSize = len(MomGrid)
TauGridSize = len(TauGrid)

shape = (Para.Order+1, MomGridSize, TauGridSize)
Data = [data.reshape(shape) for data in Data]

fig, ax = plt.subplots()

if(XType == "Mom"):
    # ??? Order 1 sigma is a delta function of tau
    oo = 2
    Data = [np.sum(d[1:oo, ...], axis=0) for d in Data]
    # y, err = Estimate(Data, Norm, lambda d: np.average(d[1, :, :], axis=1))
    y, err = Estimate(Data, Norm, lambda d: np.average(d[:, :], axis=1))
    Errorbar(MomGrid/Para.kF, y, err, fmt='o-', color="r", label=("Order-{0}".format(oo-1)) )
    ax.set_xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    ax.set_xlabel("$K$", size=size)

    x = MomGrid
    l = np.sqrt(Para.Mass2+Para.Lambda)
    # print(Para.Mass2, Para.Lambda)
    kF = Para.kF
    # Analytic Fock energy
    y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
                      (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))

    x = kF
    Mu = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
                       (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    print("Mu: ", Mu)
    # for i in range(MomGridSize):
    #     print(f"{MomGrid[i]/Para.kF:12.6f}{(y[i]-Mu)/Para.EF:12.6f}")
    ax.plot(MomGrid/Para.kF, y, "ko-", label="Fock Analytic")
    plt.title("Fock energy")

elif(XType == "Z"):

    phyFreq = np.array([-1, 1])*np.pi / Para.Beta
    Fourier = fourier.fourier(TauGrid, phyFreq, Para.Beta)

    arr = np.amin(abs(MomGrid-Para.kF))
    kFidx = np.where(abs(arr - abs(MomGrid-Para.kF)) < 1.0e-20)[0][0]

    for o in Order:
        SigmaW, Err = Estimate(Data, Norm, lambda d: Fourier.naiveT2W(
            np.sum(d[2:o+1, :, :], axis=0)))
        # print SigmaW.shape
        # print SigmaW[:, 1]-SigmaW[:, 0]
        Errorbar(MomGrid/Para.kF, 1.0-(SigmaW[:, 1].imag-SigmaW[:, 0].imag)/(2.0*np.pi/Para.Beta),
                 color=ColorList[o], label="Order {0}".format(o))
        plt.axvline(x=1.0, linestyle='--')
    ax.set_xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    ax.set_xlabel("$Ext K$", size=size)

elif(XType == "Tau"):
    arr = np.amin(abs(MomGrid-Para.kF))
    kFidx = np.where(abs(arr - abs(MomGrid-Para.kF)) < 1.0e-20)[0][0]

    klist = [kFidx]
    # N = 8
    # o = 1
    oo = 2
    Data = [np.sum(d[1:oo, ...], axis=0) for d in Data]
    for i in klist:
        q = int(i)
        Avg, Err = Estimate(Data, Norm)
        ax.errorbar(TauGrid/Para.Beta, Avg[q], yerr=Err[q], fmt='o-',
                    capthick=1, capsize=4, color=ColorList[i], label="$k={:.2f}k_F$".format(MomGrid[q]/Para.kF))

    ax.set_xlabel("$\\tau/\\beta$")
    ax.set_xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])
    ax.set_ylim([0, 1])


elif(XType == "Freq"):
    N = 5
    # o = 2
    oo = 4
    Data = [np.sum(d[1:oo, ...], axis=0) for d in Data]

    MaxFreq = 50
    phyFreq = [(freq*2.0+1.0)*np.pi /
               Para.Beta for freq in range(-MaxFreq, MaxFreq)]
    Fourier = fourier.fourier(TauGrid, phyFreq, Para.Beta)

    for i in range(N):
        q = int(i*MomGridSize/N)
        dataW = [Fourier.naiveT2W(d[q]) for d in Data]
        SigmaW, Err = Estimate(dataW, Norm)
        # ax.errorbar(phyFreq, SigmaW.real, fmt='s-',
        #             capthick=1, capsize=2, color=ColorList[2*i+1], label="$k={0}k_F$".format(MomGrid[q]/Para.kF))
        ax.errorbar(phyFreq, SigmaW.imag, fmt='o-',
                    capthick=1, capsize=2, markersize=2, label="$k={:.2f}k_F$".format(MomGrid[q]/Para.kF))
    # ax.set_xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
# plt.tight_layout()

_ = InteractiveLegend(ax)

# plt.savefig("spin_rs1_lambda1.pdf")
plt.grid()
plt.show()