#!/usr/bin/python
from utility.IO import *
import utility.fourier as fourier

# XType = "Tau"
XType = "Mom"
OrderByOrder = False
# 0: I, 1: T, 2: U, 3: S

Para = param()
Order = range(0, Para.Order+1)

Data, Norm, Step, Grids = LoadFile("./Data", "polar_pid[0-9]+.dat")

TauGrid = Grids["TauGrid"]
MomGrid = Grids["KGrid"]
TauGridSize = len(TauGrid)
MomGridSize = len(MomGrid)

shape = (Para.Order+1, MomGridSize, TauGridSize)
Data = [data.reshape(shape) for data in Data]
print Data[0].shape
Avg, Err = Estimate(Data, Norm)

# fig, ax = plt.subplots()
plt.figure()

if(XType == "Mom"):

    phyFreq = [0.0, ]
    Fourier = fourier.fourier(TauGrid, [phyFreq, ], Para.Beta)
    for o in Order:
        yList = [Fourier.naiveT2W(d[o, :, :]) for d in Data]
        y, err = Estimate(yList, Norm)
        # err = np.average(Err[o, :, :], axis=1)
        plt.errorbar(MomGrid/Para.kF, y[:, 0], yerr=err, fmt='o-', capthick=1, capsize=4,
                     color=ColorList[o], label="Order {0}".format(o))

    # x = ExtMomBin*kF
    # l = Mass2+Lambda
    # y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
    #                   (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    # ErrorPlot(ax, ExtMomBin, y, "k", ".", "Analytic")

    plt.xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    plt.xlabel("$Ext K$", size=size)

elif(XType == "Tau"):
    N = 8
    o = 1
    for i in range(N):
        q = i*MomGridSize/N
        Avg, Err = Estimate(Data, Norm)
        plt.errorbar(TauGrid/Para.Beta, Avg[o, q, :], yerr=Err[o, q, :], fmt='o-',
                     capthick=1, capsize=4, color=ColorList[i], label="$k={0}k_F$".format(MomGrid[q]/Para.kF))

    plt.xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
# plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
