import numpy as np
from scipy import integrate
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
Order = [0, 1, 2, 3, 4, 5]
SpinIndex = 2

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
EF = kF**2
Beta /= EF

Data = {}  # key: (order)
TauBin = None
ExtMomBin = None
TauBinSize = None
ExtMomBinSize = None


for order in Order:
    files = os.listdir(folder)
    Num = 0
    Norm = 0
    Data0 = None
    # if(order == 0):
    #     FileName = "vertex{0}_pid[0-9]+.dat".format(chan)
    # else:
    #     FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)

    FileName = "polar{0}_pid[0-9]+.dat".format(order)

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
                if TauBin is None:
                    TauBin = np.fromstring(line3.split(":")[1], sep=' ')
                    TauBinSize = len(TauBin)
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
    Data0 = Data0.reshape((TauBinSize, ExtMomBinSize))

    # average the angle distribution
    Data[order] = Data0

# print DataEqT[(1)]


# def Fourier(f, t):


def ErrorPlot(p, x, d, color, marker, label=None, size=4, shift=False):
    p.plot(x, d, marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)


w = 1-0.429

fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)
ColorList = ['k', 'r', 'b', 'g', 'm', 'c', 'navy', 'y']
ColorList = ColorList*40


if(XType == "Mom"):
    for o in Order:
        ErrorPlot(ax, ExtMomBin, sum(Data[o])*Beta/TauBinSize,
                  ColorList[o], 's', "Order {0}".format(o))
    ax.set_xlim([ExtMomBin[0], ExtMomBin[-1]])
    ax.set_xlabel("$Ext K$", size=size)

    # x = ExtMomBin*kF
    # l = Mass2+Lambda
    # y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((x-kF)/l)-l/kF*np.arctan((x+kF)/l) -
    #                   (l*l-x*x+kF*kF)/4.0/x/kF*np.log((l*l+(x-kF)**2)/(l*l+(x+kF)**2)))
    # ErrorPlot(ax, ExtMomBin, y, "k", ".", "Analytic")

elif(XType == "Tau"):
    N = 8
    o = 2
    for i in range(N):
        q = i*ExtMomBinSize/N
        ErrorPlot(ax, TauBin/Beta, Data[o][:, q], ColorList[i],
                  's', "k={0}".format(ExtMomBin[q]))
        ax.set_xlim([TauBin[0]/Beta, TauBin[-1]/Beta])

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
