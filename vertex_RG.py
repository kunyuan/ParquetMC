#!/usr/bin/python
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
# XType = "Mom"
XType = "Angle"
OrderByOrder = False
# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [1, 2, 3]
SpinIndex = 2
IsFullVer4 = True
IsIrreducible = False
# IsFullVer4 = False

MaxOrder = None
rs = None
Lambda = None
Mass2 = None
Beta = None
Charge2 = None
TotalStep = None
BetaStr = None
rsStr = None
ChargeStr = None
LambdaStr = None
Mass2Str = None

with open("inlist", "r") as file:
    line = file.readline()
    para = line.split(" ")
    MaxOrder = int(para[0])
    BetaStr = para[1]
    Beta = float(BetaStr)
    rsStr = para[2]
    rs = float(rsStr)
    Mass2Str = para[3]
    Mass2 = float(Mass2Str)
    LambdaStr = para[4]
    Lambda = float(LambdaStr)
    ChargeStr = para[5]
    Charge2 = float(ChargeStr)
    TotalStep = float(para[6])

Order = range(0, MaxOrder+1)

folder = "./Beta{0}_rs{1}_lambda{2}/".format(int(Beta), rs, Mass2)
# folder = "./3_Beta{0}_lambda{2}/".format(Beta, rs, Lambda)

##############   2D    ##################################
###### Bare Green's function    #########################
# kF = np.sqrt(2.0)/rs  # 2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
# Bubble = 0.15916/2  # 2D, Beta=10, rs=1
# Bubble = 0.0795775  # 2D, Beta=20, rs=1

#############  3D  ######################################
kF = (9.0*np.pi/4.0)**(1.0/3.0)/rs
Nf = kF/2.0/np.pi**2
Bubble = 0.0971916  # 3D, Beta=10, rs=1

Data = {}  # key: (order, channel)
DataWithAngle = {}  # key: (order, channel)

Para = param()
AngleBin = BuildAngleGrid(AngGridSize)
ExtMomBin = BuildMomGrid(Para.MaxExtMom, MomGridSize)
AngleBinSize = AngGridSize
ExtMomBinSize = MomGridSize


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


def Mirror(x, y):
    x2 = np.zeros(len(x)*2)
    x2[:len(x)] = x
    x2[len(x):] = -x[::-1]
    y2 = np.zeros(len(y)*2)
    y2[:len(y)] = y
    y2[len(y):] = y[::-1]
    # return x2, y2
    return np.copy(x), np.copy(y)

# def TauIntegration(Data):
#     return np.sum(Data, axis=-1) * \
#         Beta/kF**2/TauBinSize


for order in Order:
    for chan in Channel:

        files = os.listdir(folder)
        Num = 0
        Norm = 0
        Data0 = None
        # if(order == 0):
        #     FileName = "vertex{0}_pid[0-9]+.dat".format(chan)
        # else:
        #     FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)

        FileName = "vertex_pid[0-9]+.dat"

        for f in files:
            if re.match(FileName, f):
                print "Loading ", f
                with open(folder+f, "r") as file:
                    line0 = file.readline()
                    Step = int(line0.split(":")[-1])/1000000
                    # print "Step:", Step
                    line1 = file.readline()
                    # print line1
                    Norm += float(line1.split(":")[-1])
                    line3 = file.readline()
                    if AngleBin is None:
                        AngleBin = np.fromstring(line3.split(":")[1], sep=' ')
                        AngleBinSize = len(AngleBin)
                    line4 = file.readline()
                    if ExtMomBin is None:
                        ExtMomBin = np.fromstring(line4.split(":")[1], sep=' ')
                        ExtMomBinSize = len(ExtMomBin)
                        ExtMomBin /= kF
                Num += 1
                d = np.loadtxt(folder+f)
                if Data0 is None:
                    Data0 = d
                else:
                    Data0 += d
        Data0 /= Norm
        Data0 = Data0.reshape(
            (Para.Order+1, 4, AngleBinSize, ExtMomBinSize, 2))

        DataWithAngle[(order, chan)] = Data0[order, chan, ...]

        # average the angle distribution
        Data[(order, chan)] = AngleIntegation(Data0[order, chan, ...], 0)


# def ErrorPlot(p, x, d, color, marker, label=None, size=4, shift=False):
#     p.plot(x, d, marker=marker, c=color, label=label,
#            lw=1, markeredgecolor="None", linestyle="--", markersize=size)


w = 1-0.429

# fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)
# ColorList = ['k', 'r', 'b', 'g', 'm', 'c', 'navy', 'y']
# ColorList = ColorList*40

if(XType == "Scale"):
    for i in range(ExtMomBinSize/4):
        index = 4*i
        ErrorPlot(ax, ScaleBin[:-2], diffData[:-2, index],
                  ColorList[i], 's', "Q {0}".format(ExtMomBin[index]))
    ax.set_xlim([0.0, ScaleBin[-2]])
    ax.set_xlabel("$Scale$", size=size)
elif (XType == "Mom"):
    i = 0
    for chan in Channel:
        if(chan == 1):
            qData = 8.0*np.pi*Charge2/(ExtMomBin**2*kF**2+Lambda)
            # qData *= 0.0
        for order in Order[1:]:
            i += 1
            if(chan == 1):
                qData -= Data[(order, chan)]
            else:
                qData = Data[(order, chan)]

            # qData = np.sum(qData, axis=1)*Beta/kF**2/TauBinSize
            # qData0 = 8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData0
            # qData=8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData
            # print qData.shape, len(ExtMomBin)
            # print qData
            ErrorPlot(ax, ExtMomBin, qData,
                      ColorList[i], 's', "Loop {0}, Chan {1}".format(order, ChanName[chan]))

    for chan in Channel:
        if(chan == 1):
            qData = 8.0*np.pi*Charge2/(ExtMomBin**2*kF**2+Lambda)
            # qData *= 0.0
            qData -= Data[(0, chan)]
        else:
            qData = Data[(0, chan)]

        # qData = np.sum(qData, axis=1)*Beta/kF**2/TauBinSize
        # qData0 = 8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData0
        # qData=8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData
        # ErrorPlot(ax, ExtMomBin, qData,
        #           ColorList[0], 'o', "Chan {1}".format(0, ChanName[chan]))

    x = np.arange(0, 3.0, 0.001)
    y = x*0.0+Bubble
    for i in range(len(x)):
        if x[i] > 2.0:
            y[i] = Bubble*(1-np.sqrt(1-4/x[i]**2))
    y0 = 8.0*np.pi*Charge2/(x*x*kF*kF+Lambda)
    # ym=y0-y0*y0*y
    yphy = 8.0*np.pi/(x*x*kF*kF+Lambda+y*8.0*np.pi)

    # ax.plot(x, yphy, 'k-', lw=2, label="physical")
    ax.plot(x, y0, 'k-', lw=2, label="original")

    # ax.plot(x, y0*y0*y, 'r-', lw=2, label="wrong")

    ax.set_xlim([0.0, ExtMomBin[-1]])
    # ax.set_xscale("log")
    ax.set_xlabel("$q/k_F$", size=size)

elif(XType == "Angle"):

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    AngTotal = None
    for chan in Channel[0:]:
        AngData = -DataWithAngle[(0, chan)]*Nf
        if AngTotal is None:
            AngTotal = AngData
        else:
            AngTotal += AngData

        # AngHalf = np.arccos(AngleBin)/2.0
        # AngTotal[:, 0] += 8.0*np.pi/Mass2*Nf
        # AngTotal[:, 0] += -8.0*np.pi/((2.0*kF*np.sin(AngHalf))**2+Mass2)*Nf

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
