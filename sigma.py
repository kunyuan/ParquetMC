#!/usr/bin/python
from utility.IO import *
import utility.fourier as fourier

# XType = "Tau"
XType = "Mom"
# XType = "Z"
# XType = "Freq"
OrderByOrder = False
# 0: I, 1: T, 2: U, 3: S

Para = param()
Order = range(0, Para.Order+1)

Data, Norm, Step, Grid = LoadFile("./Data", "sigma_pid[0-9]+.dat")

MomGrid = Grid["KGrid"]
TauGrid = Grid["TauGrid"]
MomGridSize = len(MomGrid)
TauGridSize = len(TauGrid)

shape = (Para.Order+1, MomGridSize, TauGridSize)
Data = [data.reshape(shape) for data in Data]

fig, ax = plt.subplots()

if(XType == "Mom"):
    # Order 1 sigma is a delta function of tau
    o = 1
    yList = [np.average(d[o, :, :], axis=1) for d in Data]
    y, err = Estimate(yList, Norm)
    plt.errorbar(MomGrid, y, yerr=err, fmt='o-', capthick=1,
                 capsize=4, color=ColorList[o], label="Order {0}".format(o))
    ax.set_xlim([MomGrid[0], MomGrid[-1]])
    ax.set_xlabel("$Ext K$", size=size)

    x = MomGrid
    l = Para.Mass2+Para.Lambda
    kF = Para.kF
    y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
                      (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    ErrorPlot(ax, MomGrid, y, "k", ".", "Analytic")

elif(XType == "Z"):
    # for o in Order:
    #     ErrorPlot(ax, ExtMomBin, (DataW1[o]-DataW2[o])/(2.0*np.pi/Beta),
    #               ColorList[o], 's', "Order {0}".format(o))
    ax.set_xlim([MomGrid[0], MomGrid[-1]])
    ax.set_xlabel("$Ext K$", size=size)

elif(XType == "Tau"):
    N = 8
    o = 2
    for i in range(N):
        q = i*MomGridSize/N
        Avg, Err = Estimate(Data, Norm)
        ax.errorbar(TauGrid/Para.Beta, Avg[o, q, :], yerr=Err[o, q, :], fmt='o-',
                    capthick=1, capsize=4, color=ColorList[i], label="$k={0}k_F$".format(MomGrid[q]/Para.kF))

    ax.set_xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])


elif(XType == "Freq"):
    N = 8
    o = 2
    MaxFreq = 50
    phyFreq = [(freq*2.0+1.0)*np.pi /
               Para.Beta for freq in range(-MaxFreq, MaxFreq)]
    Fourier = fourier.fourier(TauGrid, phyFreq, Para.Beta)

    for i in range(N):
        q = i*MomGridSize/N
        dataW = [Fourier.naiveT2W(d[o, q, :]) for d in Data]
        SigmaW, Err = Estimate(dataW, Norm)
        # ax.errorbar(phyFreq, SigmaW.real, fmt='s-',
        #             capthick=1, capsize=2, color=ColorList[2*i+1], label="$k={0}k_F$".format(MomGrid[q]/Para.kF))
        ax.errorbar(phyFreq, SigmaW.imag, fmt='o-',
                    capthick=1, capsize=2, markersize=2, label="$k={0}k_F$".format(MomGrid[q]/Para.kF))
    # ax.set_xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])

plt.legend(loc=4, frameon=False, fontsize=size)
# plt.title("2D density integral")
# plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
