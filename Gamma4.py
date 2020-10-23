#!/usr/bin/env python3
from utility.IO import *
from utility.plot import *
import utility.angle as angle
import numpy as np
import sys
import argparse
# import matplotlib.pyplot as plt
# import matplotlib as mat
# mat.rcParams.update({'font.size': 16})
# mat.rcParams["font.family"] = "Times New Roman"
# size = 12


PlotType = "As_Aa"
PlotType = "Gamma_4spin"

parser = argparse.ArgumentParser("Specify some parameters.")
parser.add_argument("folder")
args = parser.parse_args()

folder = args.folder
print("Folder to plot : " + folder)



Para = param(folder)
# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
ChanColor = {0: "k", 1: "r", 2: "b", 3: "g"}
# 0: total, 1: order 1, ...
Order = range(Para.Order+1)
IsIrreducible = True

shape = (Para.Order+1, 4, Para.AngGridSize, Para.MomGridSize, 2)
Data, Norm, Step, Grid = LoadFile(folder, "vertex_pid[0-9]+.dat", shape)
# Data shape : (pid numbers, order, chan, AngleGrid, KGrid)

AngGrid = Grid["AngleGrid"]
MomGrid = Grid["KGrid"]
Angle = np.arccos(AngGrid)


def PrintInfo(Channel, Data, DataErr):
    Data = -np.copy(Data)
    DataErr = np.copy(DataErr)

    Data *= Para.Nf
    DataErr *= Para.Nf
    # print Data.shape, DataErr.shape
    print("{0}     Q/kF,    Data,    Error".format(Channel))
    print("As: {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[0], DataErr[0]))
    print("Aa:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[1], DataErr[1]))


def SpinMapping(Data):
    d = np.copy(Data)
    if PlotType == "As_Aa":
        d[..., 0] += d[..., 1]/Para.Spin
        d[..., 1] /= Para.Spin
    elif PlotType == "Gamma_4spin":
        e = d[..., 0] + d[..., 1]
        d[..., 1] = d[..., 0]
        d[..., 0] = e
    return d


def Bare(angle, Lambda):
    Bare = np.zeros([Para.AngGridSize, 2])
    if IsIrreducible == False:
        Bare[:, 0] += -8.0*np.pi/(Para.Mass2+Lambda)*Para.Nf
    Bare[:, 1] += +8.0 * np.pi / \
        ((2.0*Para.kF*np.sin(angle/2.0))**2+Para.Mass2+Lambda)*Para.Nf
    Bare = SpinMapping(Bare)
    return Bare

if PlotType == "As_Aa":
    label_bare_1 = "$u_s$"
    label_bare_2 = "$u_a$"
    label_all_1 = "$A_s$"
    label_all_2 = "$A_a$"
elif PlotType == "Gamma_4spin":
    label_bare_1 = "$u_{\\uparrow\\uparrow}$"
    label_bare_2 = "$u_{\\uparrow\\downarrow}$"
    label_all_1 = "$A_{\\uparrow\\uparrow}$"
    label_all_2 = "$A_{\\uparrow\\downarrow}$"

ColorList = ["grey", "r", "b", "g", "k", "o", "y"] * 4
bare = Bare(Angle, 0.0)
bareLambda = Bare(Angle, Para.Lambda)


# ---------------------------- Converge ----------------------------------

# fig, (cnvg1, cnvg2) = plt.subplots(1, 2, sharey=True)
# plt.suptitle("Converge: $r_s={0}, \\lambda={1}$".format(Para.Rs, Para.Lambda))
# cnvg1.plot(Angle, -bare[:, 0], '-', c='y', label="bare")
# cnvg2.plot(Angle, -bare[:, 1], '-', c='y', label="bare")
# cnvg1.plot(Angle, -bareLambda[:, 0], '-', c=ColorList[0], label="bare_$\\lambda$")
# cnvg2.plot(Angle, -bareLambda[:, 1], '-', c=ColorList[0], label="bare_$\\lambda$")
# for o in range(2, Para.Order+1):
#     label_str = "order=" + str(o)
#     data_order = [np.sum(d[1:o, ...], axis=0) for d in Data]
#     data = [SpinMapping(np.sum(d[:, :, 0, :], axis=0))*Para.Nf for d in data_order]
#     avg, err = Estimate(data, Norm)
#     bareLambda = Bare(Angle, Para.Lambda)
#     cnvg1.plot(Angle, -(avg[:, 0]+bareLambda[:, 0]), c=ColorList[o], label=label_str)
#     cnvg2.plot(Angle, -(avg[:, 1]+bareLambda[:, 1]), c=ColorList[o], label=label_str)
# cnvg2.set_ylim([-3, 2])
# plt.legend(loc=1, frameon=False, fontsize=size)

# -------------------------------------------------------------------------



# ------------------------- Gamma4,T,U,S,I ------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
plt.suptitle("$r_s={0}, \\lambda={1}$".format(Para.Rs, Para.Lambda))

ax1.plot(Angle, -bare[:, 0], '-', c='y', label=label_bare_1)
ax2.plot(Angle, -bare[:, 1], '-', c='y', label=label_bare_2)
Data = [np.sum(d[1:Para.Order+1, ...], axis=0) for d in Data]



for chan in Channel:
    if PlotType == "As_Aa":
        label_channel_1 = f"${ChanName[chan]}_s$"
        label_channel_2 = f"${ChanName[chan]}_a$"
    elif PlotType == "Gamma_4spin":
        label_channel_1 = f"${ChanName[chan]}$" + "$_{\\uparrow\\uparrow}$"
        label_channel_2 = f"${ChanName[chan]}$" + "$_{\\uparrow\\downarrow}$"

    data = [SpinMapping(d[chan, :, 0, :])*Para.Nf for d in Data]
    avg, err = Estimate(data, Norm)
    ax1.errorbar(Angle, -avg[:, 0], yerr=err[:, 0], fmt='-',
                 capthick=1, capsize=4, c=ChanColor[chan], label=label_channel_1)
    ax2.errorbar(Angle, -avg[:, 1], yerr=err[:, 1], fmt='-',
                 capthick=1, capsize=4, c=ChanColor[chan], label=label_channel_2)
    print(chan, avg[:, 0])

data = [SpinMapping(np.sum(d[:, :, 0, :], axis=0))*Para.Nf for d in Data]
avg, err = Estimate(data, Norm)
ax1.errorbar(Angle, -(avg[:, 0]+bareLambda[:, 0]), yerr=err[:, 0], fmt='-',
             capthick=1, capsize=4, c='k', label=label_all_1)
ax2.errorbar(Angle, -(avg[:, 1]+bareLambda[:, 1]), yerr=err[:, 1], fmt='-',
             capthick=1, capsize=4, c='k', label=label_all_2)

# ax1.set_xlim([-1.01, 1.01])
ax1.set_xlim([0.0, 3.15])
ax1.set_ylim([-3, 2])
ax1.set_xlabel("$\\theta$", size=size)
# ax2.set_xlim([-1.01, 1.01])
ax2.set_xlim([0.0, 3.15])
ax2.set_ylim([-3, 2])
ax2.set_xlabel("$\\theta$", size=size)
ax1.legend(loc=1, frameon=False, fontsize=size)
ax2.legend(loc=1, frameon=False, fontsize=size)
# ax1.set_xticks([0, 1, 2, 3])
# ax2.set_xticks([0, 1, 2, 3])

plt.legend(loc=1, frameon=False, fontsize=size)

plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
