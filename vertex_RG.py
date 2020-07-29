#!/usr/bin/env python3
from utility.IO import *
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

Para = param()
# XType = "Tau"
# XType = "Mom"
XType = "Angle"
OrderByOrder = False
# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = range(Para.Order+1)
SpinIndex = 2
IsFullVer4 = True
IsIrreducible = False
# IsFullVer4 = False

shape = (Para.Order+1, 4, Para.AngGridSize, Para.MomGridSize, 2)

Data, Norm, Step, Grid = LoadFile("./Data", "vertex_pid[0-9]+.dat", shape)

AngGrid = Grid["AngleGrid"]
MomGrid = Grid["KGrid"]


def SpinMapping(Data):
    d = np.copy(Data)
    d[..., 0] += d[..., 1]/Para.Spin
    d[..., 1] /= Para.Spin
    return d


def AngleIntegation(Data, l):
    # l: angular momentum
    shape = Data.shape[1:]
    Result = np.zeros(shape)
    for x in range(AngleBinSize):
        # Result += Data[x, ...] * \
        #     np.cos(l*AngleBin[x])*2.0*np.pi/AngleBinSize
        Result += Data[x, ...]*2.0/AngleBinSize
    return Result/2.0
    # return Result


fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

AngTotal = None
for chan in Channel[0:]:
    AngData = -DataWithAngle[(0, chan)]*Nf
    if AngTotal is None:
        AngTotal = AngData
    else:
        AngTotal += AngData

    x = np.arccos(AngleBin)
    y = AngData[:, 0, 0]+AngData[:, 0, 1]/SpinIndex

    print "Chan {0} As: {1}".format(ChanName[chan], np.mean(y))

    ErrorPlot(ax1, x, y, ColorList[chan], 's',
              "q/kF={0}, {1}, As, ".format(ExtMomBin[0], ChanName[chan]))

    x = np.arccos(AngleBin)
    y = AngData[:, 0, 1]/SpinIndex

    ErrorPlot(ax2, x, y, ColorList[chan], 's',
              "q/kF={0}, {1}, Aa, ".format(ExtMomBin[0], ChanName[chan]))

AngHalf = np.arccos(AngleBin)/2.0
if IsFullVer4:
    Bare = np.zeros_like(AngTotal[:, 0, :])
    if IsIrreducible == False:
        Bare[:, 0] += 8.0*np.pi/(Mass2+Lambda)*Nf
    Bare[:, 1] += -8.0 * np.pi / \
        ((2.0*kF*np.sin(AngHalf))**2+Mass2+Lambda)*Nf

    AngTotal[:, 0, 0] += Bare[:, 0]
    AngTotal[:, 0, 1] += Bare[:, 1]

    x = np.arccos(AngleBin)
    ErrorPlot(ax1, x, Bare[:, 0]+Bare[:, 1]/SpinIndex, ColorList[-1], 's',
              "q/kF={0}, Bare".format(ExtMomBin[0]))

    x = np.arccos(AngleBin)
    ErrorPlot(ax2, x, Bare[:, 1]/SpinIndex, ColorList[-1], 's',
              "q/kF={0}, Bare".format(ExtMomBin[0]))

AngLandau = AngTotal+DataWithAngle[(0, 1)]*Nf

if SpinIndex == 2:
    x = np.arccos(AngleBin)
    y = AngTotal[:, 0, 0]+AngTotal[:, 0, 1]/2.0
    ErrorPlot(ax3, x, y, ColorList[0], 's',
              "q/kF={0}, As".format(ExtMomBin[0]))
    print "As: ", sum(y)/len(y)

    x = np.arccos(AngleBin)
    y = AngTotal[:, 0, 1]/2.0
    ErrorPlot(ax3, x, y, ColorList[1], 's',
              "q/kF={0}, Aa".format(ExtMomBin[0]))
    print "Aa: ", sum(y)/len(y)

    x = np.arccos(AngleBin)
    y = AngLandau[:, 0, 0]+AngLandau[:, 0, 1]/2.0
    ErrorPlot(ax3, x, y, ColorList[2], 's',
              "q/kF={0}, Fs".format(ExtMomBin[0]))
    print "Fs: ", sum(y)/len(y)

    x = np.arccos(AngleBin)
    y = AngLandau[:, 0, 1]/2.0
    ErrorPlot(ax3, x, y, ColorList[3], 's',
              "q/kF={0}, Fa".format(ExtMomBin[0]))
    print "Fa: ", sum(y)/len(y)
else:
    x = np.arccos(AngleBin)
    y = AngTotal[:, 0, 0]+AngTotal[:, 0, 1]
    ErrorPlot(ax3, x, y, ColorList[-1], 's',
              "q/kF={0}, A".format(ExtMomBin[0]))

# ax.set_xlim([-np.arccos(AngleBin[0]), np.arccos(AngleBin[0])])
ax1.set_xlim([0.0, np.pi])
# ax.set_ylim([0.0, 5.0])
ax1.set_xlabel("$Angle$", size=size)
ax2.set_xlim([0.0, np.pi])
ax3.set_xlim([0.0, np.pi])
ax2.set_xlabel("$Angle$", size=size)
ax3.set_xlabel("$Angle$", size=size)
ax1.legend(loc=1, frameon=False, fontsize=size)
ax2.legend(loc=1, frameon=False, fontsize=size)
ax3.legend(loc=1, frameon=False, fontsize=size)
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([-0.02, 0.125])
# ax.set_ylim([0.07, 0.125])
# ax.xaxis.set_label_coords(0.97, -0.01)
# # ax.yaxis.set_label_coords(0.97, -0.01)
# ax.text(-0.012,0.52, "$-I$", fontsize=size)
# ax.set_ylabel("$-\Gamma_4(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

# plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
# plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
