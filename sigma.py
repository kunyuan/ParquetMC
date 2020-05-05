import numpy as np
from scipy import integrate
from utility import *
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob
import os
import re
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

XType = "Tau"
# XType = "Mom"
# XType = "Z"
OrderByOrder = False
SpinIndex = 2
Dim = 3
# 0: I, 1: T, 2: U, 3: S

Para = param("inlist", Dim, SpinIndex)
Order = range(0, Para.Order+1)
TauBin, ExtMomBin, TauBinSize, ExtMomBinSize = [None, ]*4
Data = {}  # key: (order)

folder = "./Beta{0}_rs{1}_lambda{2}/".format(
    int(Para.Beta*Para.EF), Para.Rs, Para.Mass2)
files = getListOfFiles(folder)

for order in Order:
    Norm = 0.0
    Data[order] = None
    for f in files:
        if re.search("sigma{0}_pid[0-9]+.dat".format(order), f):
            print "Loading ", f
            with open(f, "r") as file:
                Step = int(file.readline().split(":")[-1])/1000000
                Norm += float(file.readline().split(":")[-1])
                TauBin = np.fromstring(file.readline().split(":")[1], sep=' ')
                TauBinSize = len(TauBin)
                ExtMomBin = np.fromstring(
                    file.readline().split(":")[1], sep=' ')
                ExtMomBinSize = len(ExtMomBin)
            d = np.loadtxt(f)
            if Data[order] is None:
                Data[order] = d
            else:
                Data[order] += d
    Data[order] /= Norm
    Data[order] = Data[order].reshape((TauBinSize, ExtMomBinSize))

fig, ax = plt.subplots()

if(XType == "Mom"):
    # Order 1 sigma is a delta function of tau
    ErrorPlot(ax, ExtMomBin, Data[o][0, :]*Para.Beta/TauBinSize,
              ColorList[o], 's', "Order {0}".format(o))
    ax.set_xlim([ExtMomBin[0], ExtMomBin[-1]])
    ax.set_xlabel("$Ext K$", size=size)

    x = ExtMomBin
    l = Para.Mass2+Para.Lambda
    kF = Para.kF
    y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
                      (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    ErrorPlot(ax, ExtMomBin, y, "k", ".", "Analytic")

elif(XType == "Z"):
    # for o in Order:
    #     ErrorPlot(ax, ExtMomBin, (DataW1[o]-DataW2[o])/(2.0*np.pi/Beta),
    #               ColorList[o], 's', "Order {0}".format(o))
    ax.set_xlim([ExtMomBin[0], ExtMomBin[-1]])
    ax.set_xlabel("$Ext K$", size=size)

elif(XType == "Tau"):
    N = 8
    o = 2
    for i in range(N):
        q = i*ExtMomBinSize/N
        ErrorPlot(ax, TauBin/Para.Beta, Data[o][:, q], ColorList[i],
                  's', "k={0}".format(ExtMomBin[q]))
        ax.set_xlim([TauBin[0]/Para.Beta, TauBin[-1]/Para.Beta])

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
