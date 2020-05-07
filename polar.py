import numpy as np
from scipy import integrate
from utility import *
from grid import *
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob
import os
import re
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

# XType = "Tau"
XType = "Mom"
OrderByOrder = False
# 0: I, 1: T, 2: U, 3: S

Para = param()
Order = range(0, Para.Order+1)
TauGrid = BuildTauGrid(Para.Beta, TauGridSize)
MomGrid = BuildMomGrid(Para.MaxExtMom, MomGridSize)

folder = "./Beta{0}_rs{1}_lambda{2}/".format(
    int(Para.Beta*Para.EF), Para.Rs, Para.Mass2)

Data = {}

for order in Order:
    Norm = 0.0
    Data[order] = None
    for f in files:
        if re.search("polar{0}_pid[0-9]+.dat".format(order), f):
            print "Loading ", f
            with open(f, "r") as file:
                Step = int(file.readline().split(":")[-1])/1000000
                Norm += float(file.readline().split(":")[-1])
            d = np.loadtxt(f)
            if Data[order] is None:
                Data[order] = d
            else:
                Data[order] += d
    Data[order] /= Norm
    Data[order] = Data[order].reshape((TauBinSize, ExtMomBinSize))


fig, ax = plt.subplots()

if(XType == "Mom"):
    for o in Order:
        print Data[o].shape
        ErrorPlot(ax, MomGrid, np.sum(Data[o][:, :], axis=0)*Para.Beta/TauGridSize,
                  ColorList[o], 's', "Order {0}".format(o))
    ax.set_xlim([MomGrid[0], MomGrid[-1]])
    ax.set_xlabel("$Ext K$", size=size)

    # x = ExtMomBin*kF
    # l = Mass2+Lambda
    # y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
    #                   (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    # ErrorPlot(ax, ExtMomBin, y, "k", ".", "Analytic")

elif(XType == "Tau"):
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
