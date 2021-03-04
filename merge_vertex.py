#!/usr/bin/env python3
from utility.IO import *
import time
import traceback
import numpy as np
import utility.angle as legendre
import sys
import argparse
import KOinteraction as KO

parser = argparse.ArgumentParser("Specify some parameters.")
parser.add_argument("folder")
args = parser.parse_args()

folder = args.folder
print("Folder to plot : " + folder)


save = True
Type = "Symmetry_Anti"
# Type = "Gamma_4spin"

SleepTime = 5
Irreducible = False

Para = param(folder)

# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = range(Para.Order+1)
# Order = [0, 1, 2, 3, ]

shape = (Para.Order+1, 4, Para.AngGridSize, Para.MomGridSize, 2)


def SpinMapping(Data):
    d = np.copy(Data)
    if Type == "Symmetry_Anti":
        d[..., 0] += d[..., 1]/Para.Spin
        d[..., 1] /= Para.Spin
    elif Type == "Gamma_4spin":
        e = d[..., 0] + d[..., 1]
        d[..., 1] = d[..., 0]
        d[..., 0] = e
    return d


def PrintInfo(Channel, Data, DataErr):
    Data = -np.copy(Data)
    DataErr = np.copy(DataErr)

    coeff = GetCoeff()
    Data *= coeff
    DataErr *= coeff

    # print Data.shape, DataErr.shape

    print("{0}     Q/kF,    Data,    Error".format(Channel))
    print("As: {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[0], DataErr[0]))
    print("Aa:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[1], DataErr[1]))


def GetCoeff():
    m = 1.0/2.0
    Zname = "ZandM_order{0}.data".format(Para.Order-1)
    try:
        Z, mStar = np.loadtxt(os.path.join(folder, "selfconsistent", Zname))
    except Exception as e:
        print("Can not load Z and m*, will assume Z=1 and m=1/2")
        Z, mStar = 1, m
    # Z, mStar = 0.873, 0.955*0.5
    coeff = Z*Z * mStar * Para.kF / (np.pi*np.pi)
    return coeff


# rs=1  Z=0.873, m*=0.95
while True:

    time.sleep(SleepTime)

    try:
        Data, Norm, Step, Grid = LoadFile(
            folder, "vertex_pid[0-9]+.dat", shape)

        AngGrid = Grid["AngleGrid"]
        MomGrid = Grid["KGrid"]

        DataAngle, ErrAngle = Estimate(Data, Norm)
        print("Write Weight file.")
        with open("weight.data", "w") as file:
            for chan in Channel:
                for angle in range(Para.AngGridSize):
                    for qidx in range(Para.MomGridSize):
                        for Dir in range(2):
                            file.write("{0} ".format(
                                DataAngle[0, chan, angle, qidx, Dir]))

        # with open("data.data", "a") as file:
        #     file.write("Dir: {0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(
        #         Data[(0, 1)][0, 0], Data[(0, 1)][0, 0], Data[(0, 2)][0, 0], Data[(0, 3)][0, 0]))
        #     file.write("Ex: {0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(
        #         Data[(0, 1)][0, 1], Data[(0, 1)][0, 1], Data[(0, 2)][0, 1], Data[(0, 3)][0, 1]))

        # Keep the ExtMom=0 elements only, and average the angle !!!!!!!
        # DataList = [np.average(d[:, :, :, 0, :], axis=2) for d in Data]

        DataList = [legendre.LegendreCoeff(d[:, :, :, 0, :], AngGrid, [
                                           0, ], axis=2)[0] for d in Data]

        # DataList = [d[:, :, 0, 0, :] for d in DataList]

        # construct bare interaction
        Bare = np.zeros(2)
        # if not Irreducible:
        #     Bare[0] -= 8.0*np.pi/(Para.Mass2+Para.Lambda)

        AngHalf = np.arccos(AngGrid)/2.0
        # ExBare = +8.0 * np.pi / \
        #     ((2.0*Para.kF*np.sin(AngHalf))**2+Para.Mass2+Para.Lambda)
        ExRs, ExRa = KO.InterFreq(
            np.array([0]), 2.0*Para.kF*np.sin(AngHalf), Para, addBare=True)
        ExR = ExRs-ExRa

        # ExBare = +8.0 * np.pi / \
        # ((2.0*Para.kF*np.sin(AngHalf))**2+Para.Mass2)
        # print ExBare.shape

        # print "ExBare: ", AngleIntegation(ExBare, 0)
        # Bare[1] = np.average(ExBare)
        Bare[1] = legendre.LegendreCoeff(ExR, AngGrid, [0, ], 0)[0]
        ExRs = legendre.LegendreCoeff(ExRs, AngGrid, [0, ], 0)[0]
        ExRa = legendre.LegendreCoeff(ExRa, AngGrid, [0, ], 0)[0]
        # ExRs *= Para.Nf
        # ExRa *= Para.Nf
        Bare[0] = ExRa*2
        Bare[1] = ExRs-ExRa
        # print(ExRs*Para.Nf, ExRa*Para.Nf)
        exchange0 = 8.0*np.pi/2.0/Para.kF**2 * \
            np.log((Para.Mass2+Para.Lambda+4.0*Para.kF**2) /
                   (Para.Mass2+Para.Lambda))/2.0  # factor 2 comes from the normalization
        print(f"Benchmark exchange bare for l=0: {Bare[1]} vs {exchange0}")
        Bare = SpinMapping(Bare)
        # print(Bare*Para.Nf)

        # print(AngGrid)
        # print(ExBare)
        # print(np.sum(ExBare)/len(ExBare), np.average(ExBare))
        # coeff = legendre.LegendreCoeff(ExBare, AngGrid, [0, ], 0)
        # print(f"{Bare[1]} vs {coeff}")

        # Bare *= 0.0
        # print Bare
        dataStrAs = ""
        dataStrAa = ""
        for o in Order[1:]:
            print(green("Order {0}".format(o)))

            # sum all orders
            # DataAllList = [np.sum(d[1:o+1, ...], axis=0) for d in DataList]
            DataAllList = [np.sum(d[o:o+1, ...], axis=0) for d in DataList]
            # sum all four channels
            DataAllList = [np.sum(d, axis=0) for d in DataAllList]
            # map DIR, EX to As, Aa
            DataAllList = [SpinMapping(d) for d in DataAllList]
            Data, Err = Estimate(DataAllList, Norm)
            # Data += Bare  # I channel has a bare part
            PrintInfo("Sum", Data, Err)

            # # print the S channel:
            # DataAllList = [np.sum(d[1:o+1, ...], axis=0) for d in DataList]
            # # sum all four channels
            # DataAllList = [d[3, ...] for d in DataAllList]
            # # map DIR, EX to As, Aa
            # DataAllList = [SpinMapping(d) for d in DataAllList]
            # Data, Err = Estimate(DataAllList, Norm)
            # PrintInfo("S channel: ", Data, Err)

            if save and o == Para.Order:
                coeff = GetCoeff()
                DataSave = Data * coeff
                ErrSave = Err * coeff
                dataStrAs += "{0:10.6f}  {1:10.6f}  ".format(
                    DataSave[0], ErrSave[0])
                dataStrAa += "{0:10.6f}  {1:10.6f}  ".format(
                    DataSave[1], ErrSave[1])
            # fig, ax1 = plt.subplots(1)
            # ax1.errorbar(MomGrid, Data[], yerr=err[:, 0], fmt='-',
            #     capthick=1, capsize=4, c=ChanColor[chan], label=f"${ChanName[chan]}_s$")

        #     # qData = Data[(o, 1)]
        #     # qDataErr = DataErr[(o, 1)]
        #     # PrintInfo("I", Data[(o, 0)], DataErr[(o, 0)])
        #     # PrintInfo("T", Data[(o, 1)], DataErr[(o, 1)])
        #     # PrintInfo("U", Data[(o, 2)], DataErr[(o, 2)])
        #     # PrintInfo("S", Data[(o, 3)], DataErr[(o, 3)])
        #     PrintInfo("Sum", qData, qDataErr)
        #     # print "\n"

        if save:
            with open("vertexq0.data", "a") as f:
                f.write("{0:4.2f} {1:4.2f} {2:4.2f} {3:4.2f} {4:4.2f} ".format(
                    Para.Beta*Para.EF, Para.Rs, Para.Mass2, Para.Lambda, Para.Order))
                f.write(dataStrAs + dataStrAa + "\n")

        print("\n")
        flag = np.array([step/1000000 >= Para.TotalStep for step in Step])
        if np.all(flag == True):
            print("End of Simulation!")
            sys.exit(0)
        sys.exit(0)

    except Exception as e:
        print(e)
        traceback.print_exc()
        sys.exit(0)
        pass
