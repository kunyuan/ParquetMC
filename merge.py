import time
from color import *
from utility import *
from grid import *

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


def AngleIntegation(Data, l):
    # l: angular momentum
    shape = Data.shape[1:]
    Result = np.zeros(shape)
    for x in range(AngGrid):
        # Result += Data[x, ...] * \
        #     np.cos(l*AngleBin[x])*2.0*np.pi/AngleBinSize
        Result += Data[x, ...]*2.0/AngGrid
    return Result/2.0
    # return Result


def SpinMapping(Data):
    d = np.copy(AngleIntegation(Data, 0))
    d[..., 0] += d[..., 1]/SpinIndex
    d[..., 1] /= SpinIndex
    return d


def PrintInfo(Channel, Data, DataErr):
    i = 0
    Data = -np.copy(Data)
    DataErr = np.copy(DataErr)

    Data[:, 0] *= Para.Nf
    Data[:, 1] *= Para.Nf
    DataErr[:, 0] *= Para.Nf
    DataErr[:, 1] *= Para.Nf

    print "{0}     Q/kF,    Data,    Error".format(Channel)
    qData0 = Data[:, 0]
    print "As: {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[i], qData0[i], DataErr[i, 0])
    qData1 = Data[:, 1]
    print "Aa:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[i], qData1[i], DataErr[i, 1])

    # print "As:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
    #     ExtMomBin[i], qData0[i]+qData1[i]/SpinIndex, DataErr[i, 0]+DataErr[i, 1]/SpinIndex)
    # print "Aa:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
    #     ExtMomBin[i], qData1[i]/SpinIndex, DataErr[i, 1]/SpinIndex)


while True:

    time.sleep(SleepTime)

    DataAngle, ErrAngle, Step = LoadFile(folder, filename, shape)
    Data, Err, Step = LoadFile(folder, filename, shape, SpinMapping)

    if len(Data) > 0:
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

        # construct bare interaction
        AngHalf = np.arccos(AngGrid)/2.0
        Bare = np.zeros(2)
        if IsIrreducible == False:
            Bare[0] -= 8.0*np.pi/(Para.Mass2+Para.Lambda)
        ExBare = +8.0 * np.pi / \
            ((2.0*Para.kF*np.sin(AngHalf))**2+Para.Mass2+Para.Lambda)
        # print "ExBare: ", AngleIntegation(ExBare, 0)
        Bare[1] += AngleIntegation(ExBare, 0)
        Bare = SpinMapping(Bare)

        qData = np.zeros_like(Data[1, 0])
        qData[0, :] += Bare[:]
        qDataErr = np.zeros_like(DataErr[(1, 0)])
        for o in Order[1:]:
            print green("Order {0}".format(o))
            qData += sum([Data[(o, i)] for i in range(4)])
            qDataErr += sum([DataErr[(o, i)] for i in range(4)])

            # qData += Data[(o, 1)]
            # qDataErr += DataErr[(o, 1)]

            # qData = Data[(o, 1)]
            # qDataErr = DataErr[(o, 1)]
            # PrintInfo("I", Data[(o, 0)], DataErr[(o, 0)])
            # PrintInfo("T", Data[(o, 1)], DataErr[(o, 1)])
            # PrintInfo("U", Data[(o, 2)], DataErr[(o, 2)])
            # PrintInfo("S", Data[(o, 3)], DataErr[(o, 3)])
            PrintInfo("Sum", qData, qDataErr)
            # print "\n"

        qData = sum([Data[(0, i)] for i in range(4)])
        qData[0, :] += Bare[:]
        qDataErr = sum([DataErr[(0, i)] for i in range(4)])
        PrintInfo("\nTotal", qData, qDataErr)

    if Step >= TotalStep:
        print "End of Simulation!"
        sys.exit(0)
