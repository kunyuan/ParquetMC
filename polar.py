from scipy import integrate
from utility import *

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
Avg, Err = Estimate(Data, Norm)

# fig, ax = plt.subplots()
plt.figure()

if(XType == "Mom"):
    for o in Order:
        yList = [np.average(d[o, :, :], axis=1) for d in Data]
        y, err = Estimate(yList, Norm)
        # err = np.average(Err[o, :, :], axis=1)
        plt.errorbar(MomGrid/Para.kF, y, yerr=err, fmt='o-', capthick=1, capsize=4,
                     color=ColorList[o], label="Order {0}".format(o))

    plt.xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    plt.xlabel("$Ext K$", size=size)

    # x = ExtMomBin*kF
    # l = Mass2+Lambda
    # y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
    #                   (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    # ErrorPlot(ax, ExtMomBin, y, "k", ".", "Analytic")

elif(XType == "Tau"):
    pass
    # N = 8
    # o = 2
    # for i in range(N):
    #     q = i*ExtMomBinSize/N
    #     ErrorPlot(ax, TauBin/Beta, Data[o][:, q], ColorList[i],
    #               's', "k={0}".format(ExtMomBin[q]))
    #     ax.set_xlim([TauBin[0]/Beta, TauBin[-1]/Beta])

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
