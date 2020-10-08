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

AngGrid = Grid["AngleGrid"]
MomGrid = Grid["KGrid"]
Angle = np.arccos(AngGrid)

# Data shape : (pid numbers, order, chan, AngleGrid, KGrid)

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
    d[..., 0] += d[..., 1]/Para.Spin
    d[..., 1] /= Para.Spin
    return d


def Bare(angle, Lambda):

    Bare = np.zeros([Para.AngGridSize, 2])
    if IsIrreducible == False:
        Bare[:, 0] += -8.0*np.pi/(Para.Mass2+Lambda)*Para.Nf
    Bare[:, 1] += +8.0 * np.pi / \
        ((2.0*Para.kF*np.sin(angle/2.0))**2+Para.Mass2+Lambda)*Para.Nf
    Bare = SpinMapping(Bare)
    return Bare

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

bare = Bare(Angle, 0.0)
ax1.plot(Angle, -bare[:, 0], '-', c='y', label="$u_s$")
ax2.plot(Angle, -bare[:, 1], '-', c='y', label="$u_a$")


Data = [np.sum(d[1:Para.Order+1, ...], axis=0) for d in Data]

# sys.exit()

# DataAllList = [np.sum(d, axis=0) for d in Data]
# map DIR, EX to As, Aa
# DataAllList = [SpinMapping(d) for d in DataAllList]
# Data, Err = Estimate(DataAllList, Norm)
# PrintInfo("Sum", Data, Err)

ColorList = ["grey", "r", "b", "g"]
for chan in Channel:
    data = [SpinMapping(d[chan, :, 0, :])*Para.Nf for d in Data]
    avg, err = Estimate(data, Norm)
    ax1.errorbar(Angle, -avg[:, 0], yerr=err[:, 0], fmt='-',
                 capthick=1, capsize=4, c=ChanColor[chan], label=f"${ChanName[chan]}_s$")
    ax2.errorbar(Angle, -avg[:, 1], yerr=err[:, 1], fmt='-',
                 capthick=1, capsize=4, c=ChanColor[chan], label=f"${ChanName[chan]}_a$")
    # print(chan, avg[:, 0])

data = [SpinMapping(np.sum(d[:, :, 0, :], axis=0))*Para.Nf for d in Data]
avg, err = Estimate(data, Norm)
bareLambda = Bare(Angle, Para.Lambda)
ax1.errorbar(Angle, -(avg[:, 0]+bareLambda[:, 0]), yerr=err[:, 0], fmt='-',
             capthick=1, capsize=4, c='k', label=f"$A_s$")
ax2.errorbar(Angle, -(avg[:, 1]+bareLambda[:, 1]), yerr=err[:, 1], fmt='-',
             capthick=1, capsize=4, c='k', label=f"$A_a$")

ax1.set_xlim([0.0, np.pi])
# ax.set_ylim([0.0, 5.0])
ax1.set_xlabel("$\\theta$", size=size)
ax2.set_xlim([0.0, np.pi])
ax2.set_xlabel("$\\theta$", size=size)
ax1.legend(loc=1, frameon=False, fontsize=size)
ax2.legend(loc=1, frameon=False, fontsize=size)
ax1.set_xticks([0, 1, 2, 3])
ax2.set_xticks([0, 1, 2, 3])

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
