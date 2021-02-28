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



# XType = "Tau"
XType = "Mom"
OrderByOrder = False
# Unit=1.0
# 0: I, 1: T, 2: U, 3: S

Para = param(folder)
Order = range(0, Para.Order+1)

Data, Norm, Step, Grids = LoadFile(folder, "polar_pid[0-9]+.dat")

TauGrid = Grids["TauGrid"]
MomGrid = Grids["KGrid"]
TauGridSize = len(TauGrid)
MomGridSize = len(MomGrid)

shape = (Para.Order+1, MomGridSize, TauGridSize)
Data = [data.reshape(shape) for data in Data]
Avg, Err = Estimate(Data, Norm)

fig, ax = plt.subplots()
# print "original ", fig
# fig, ax = plt.figure()
lines = []

if(XType == "Mom"):

    phyFreq = [0.0, ]
    Fourier = fourier.fourier(TauGrid, [phyFreq, ], Para.Beta)
    for o in Order[1:]:
        y, err = Estimate(Data, Norm, lambda x: Fourier.naiveT2W(
            np.sum(x[1:o+1, :, :], axis=0)))
        Errorbar(MomGrid/Para.kF, y[:, 0], err[:, 0]*2.0, fmt='o-',
                 color=ColorList[o], label="Order {0}".format(o))
        print("order:{0}   {1} +- {2}".format(o, y[0, 0], 2*err[0] ))
    # x = ExtMomBin*kF
    # l = Mass2+Lambda
    # y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
    #                   (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    # ErrorPlot(ExtMomBin, y, y*0.0, "k", ".", "Analytic")

    ax.set_xlim([MomGrid[0]/Para.kF, MomGrid[-1]/Para.kF])
    ax.set_xlabel("$Ext K$", size=size)

elif(XType == "Tau"):
    N = 8
    o = 1
    for i in range(N):
        q = int(i*MomGridSize/N)
        Avg, Err = Estimate(Data, Norm)
        plt.errorbar(TauGrid/Para.Beta, Avg[o, q, :], yerr=Err[o, q, :], fmt='o-',
                     capthick=1, capsize=4, color=ColorList[i], label="$k={0}k_F$".format(MomGrid[q]/Para.kF))

    plt.xlim([TauGrid[0]/Para.Beta-1e-3, TauGrid[-1]/Para.Beta])

# leg = ax.legend(handler_map=my_handler_map)
leg = ax.legend(loc=1, frameon=False, fontsize=size)

_ = InteractiveLegend(ax)  # keep the obj to prevent be collected by gc

# plt.savefig("spin_rs1_lambda1.pdf")
# fig.canvas.mpl_connect('pick_event', l.onpick)
plt.grid()
plt.show()
