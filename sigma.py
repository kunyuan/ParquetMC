from scipy import integrate
from utility import *
from grid import *

XType = "Tau"
# XType = "Mom"
# XType = "Z"
OrderByOrder = False
# 0: I, 1: T, 2: U, 3: S

Para = param()
Order = range(0, Para.Order+1)
TauGrid = BuildTauGrid(Para.Beta, TauGridSize)
MomGrid = BuildMomGrid(Para.MaxExtMom, MomGridSize)

folder = "./Beta{0}_rs{1}_lambda{2}/".format(
    int(Para.Beta*Para.EF), Para.Rs, Para.Mass2)

filename = "sigma_pid[0-9]+.dat"

shape = (Para.Order+1, MomGridSize, TauGridSize)

Avg, Err, Step = LoadFile(folder, filename, shape)

fig, ax = plt.subplots()

if(XType == "Mom"):
    # Order 1 sigma is a delta function of tau
    o = 1
    y = np.average(Avg[o, :, :], axis=1)
    err = np.average(Err[o, :, :], axis=1)/np.sqrt(len(TauGrid))
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
        ax.errorbar(TauGrid/Para.Beta, Avg[o, q, :], yerr=Err[o, q, :], fmt='o-',
                    capthick=1, capsize=4, color=ColorList[i], label="k={0}".format(MomGrid[q]))
    ax.set_xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])

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
