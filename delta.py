from scipy import integrate
from utility.IO import *

XType = "Tau"
# XType = "Mom"

Para = param()
Order = range(0, Para.Order+1)

shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)
Data, Norm, Step, Grids = LoadFile("./Data", "delta_pid[0-9]+.dat", shape)
TauGrid = Grids["TauGrid"]
MomGrid = Grids["KGrid"]

fig, ax = plt.subplots()

if(XType == "Mom"):
    # Order 1 sigma is a delta function of tau
    o = 1
    # average over the external Tau variable
    # note the 1st order is a delta(t) function, thus is indepdent of external Tau
    yList = [np.average(d[o, :, :], axis=1) for d in Data]
    y, err = Estimate(yList, Norm)
    plt.errorbar(MomGrid, y, yerr=err, fmt='o-', capthick=1,
                 capsize=4, color=ColorList[o], label="Order {0}".format(o))
    ax.set_xlim([MomGrid[0], MomGrid[-1]])
    ax.set_xlabel("$Ext K$", size=size)

    # x = MomGrid
    # l = Para.Mass2+Para.Lambda
    # kF = Para.kF
    # y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
    #                   (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    # ErrorPlot(ax, MomGrid, y, "k", ".", "Analytic")

elif(XType == "Tau"):
    N = 8
    o = 2
    for i in range(N):
        q = i*Para.MomGridSize/N
        dataList = [d[o, q, :] for d in Data]
        Avg, Err = Estimate(dataList, Norm)
        if i == N/2:
            print Avg[0], Err[0]
            for d, norm, step in zip(dataList, Norm, Step):
                print d[0]/norm, norm, step
        ax.errorbar(TauGrid/Para.Beta, Avg, yerr=Err, fmt='o-',
                    capthick=1, capsize=2, markersize=2, label="k={0}".format(MomGrid[q]/Para.kF))
    ax.set_xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])


plt.legend(loc=1, frameon=False, fontsize=size)
# plt.tight_layout()

plt.show()
