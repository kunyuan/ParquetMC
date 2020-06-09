#!/usr/bin/python
import time
from color import *
from utility import *
from grid import *
import traceback

SleepTime = 5
IsIrreducible = False

Para = param()

# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = range(Para.Order+1)
# Order = [0, 1, 2, 3, ]

folder = "./Beta{0}_rs{1}_lambda{2}/".format(
    int(Para.Beta*Para.EF), Para.Rs, Para.Mass2)

filename = "vertex_pid[0-9]+.dat"

shape = (Para.Order+1, 4, AngGridSize, MomGridSize, 2)

AngGrid = BuildAngleGrid(AngGridSize)
MomGrid = BuildMomGrid(Para.MaxExtMom, MomGridSize)


def SpinMapping(Data):
    d = np.copy(Data)
    d[..., 0] += d[..., 1]/SpinIndex
    d[..., 1] /= SpinIndex
    return d


def PrintInfo(Channel, Data, DataErr):
    Data = -np.copy(Data)
    DataErr = np.copy(DataErr)

    Data *= Para.Nf
    DataErr *= Para.Nf

    # print Data.shape, DataErr.shape

    print "{0}     Q/kF,    Data,    Error".format(Channel)
    print "As: {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[0], DataErr[0])
    print "Aa:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[1], DataErr[1])


while True:

    time.sleep(SleepTime)

    try:
        DataList, NormList, Step = LoadFile(folder, filename, shape)

        DataAngle, ErrAngle = Estimate(DataList, NormList)
        print "Write Weight file."
        with open("weight.data", "w") as file:
            for chan in Channel:
                for angle in range(AngGridSize):
                    for qidx in range(MomGridSize):
                        for Dir in range(2):
                            file.write("{0} ".format(
                                DataAngle[0, chan, angle, qidx, Dir]))

        # with open("data.data", "a") as file:
        #     file.write("Dir: {0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(
        #         Data[(0, 1)][0, 0], Data[(0, 1)][0, 0], Data[(0, 2)][0, 0], Data[(0, 3)][0, 0]))
        #     file.write("Ex: {0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(
        #         Data[(0, 1)][0, 1], Data[(0, 1)][0, 1], Data[(0, 2)][0, 1], Data[(0, 3)][0, 1]))

        # Keep the ExtMom=0 elements only, and average the angle
        DataList = [np.average(d[:, :, :, 0, :], axis=2) for d in DataList]
        # DataList = [d[:, :, 0, 0, :] for d in DataList]

        # construct bare interaction
        Bare = np.zeros(2)
        if IsIrreducible == False:
            Bare[0] -= 8.0*np.pi/(Para.Mass2+Para.Lambda)

        AngHalf = np.arccos(AngGrid)/2.0
        ExBare = +8.0 * np.pi / \
            ((2.0*Para.kF*np.sin(AngHalf))**2+Para.Mass2+Para.Lambda)
        # print ExBare.shape

        # print "ExBare: ", AngleIntegation(ExBare, 0)
        Bare[1] = np.average(ExBare)
        Bare = SpinMapping(Bare)

        # Bare *= 0.0
        # print Bare
        for o in Order[1:]:
            print green("Order {0}".format(o))

            # sum all orders
            DataAllList = [np.sum(d[1:o+1, ...], axis=0) for d in DataList]
            # sum all four channels
            DataAllList = [np.sum(d, axis=0) for d in DataAllList]
            # map DIR, EX to As, Aa
            DataAllList = [SpinMapping(d) for d in DataAllList]
            Data, Err = Estimate(DataAllList, NormList)
            Data += Bare  # I channel has a bare part
            PrintInfo("Sum", Data, Err)

        #     # qData = Data[(o, 1)]
        #     # qDataErr = DataErr[(o, 1)]
        #     # PrintInfo("I", Data[(o, 0)], DataErr[(o, 0)])
        #     # PrintInfo("T", Data[(o, 1)], DataErr[(o, 1)])
        #     # PrintInfo("U", Data[(o, 2)], DataErr[(o, 2)])
        #     # PrintInfo("S", Data[(o, 3)], DataErr[(o, 3)])
        #     PrintInfo("Sum", qData, qDataErr)
        #     # print "\n"

        print "\n"
        flag = np.array([step/1000000 >= Para.TotalStep for step in Step])
        if np.all(flag == True):
            print "End of Simulation!"
            sys.exit(0)

    except Exception as e:
        print e
        traceback.print_exc()
        pass
