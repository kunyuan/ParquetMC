import os
import sys
import re
import glob
import time
import numpy as np

SleepTime = 10
SpinIndex = 2

rs = None
Lambda = None
Mass2 = None
Beta = None
Charge2 = None
TotalStep = None
BetaStr = None
rsStr = None
ChargeStr = None
Mass2Str = None
LambdaStr = None

with open("inlist", "r") as file:
    line = file.readline()
    para = line.split(" ")
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
    TotalStep = float(para[7])

print rs, Beta, Mass2, Lambda, TotalStep

# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [0, ]

folder = "./Beta{0}_rs{1}_lambda{2}/".format(BetaStr, rsStr, Mass2Str)

AngleBin = None
ExtMomBin = None
AngleBinSize = None
ExtMomBinSize = None
Data = {}  # key: (order, channel)
DataWithAngle = {}  # key: (order, channel)
DataErr = {}  # key: (order, channel)

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
Step = None


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


def PrintInfo(Channel, Data, DataErr):
    i = 0
    Data = np.copy(Data)
    DataErr = np.copy(DataErr)
    Data[:, 0] *= Nf
    Data[:, 1] *= Nf
    DataErr[:, 0] *= Nf
    DataErr[:, 1] *= Nf
    print "\n{0}     Q/kF,    Data,    Error".format(Channel)
    qData0 = Data[:, 0]
    print "Dir: {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        ExtMomBin[i], qData0[i], DataErr[i, 0])
    qData1 = Data[:, 1]
    print "Ex:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        ExtMomBin[i], qData1[i], DataErr[i, 1])

    if SpinIndex == 2:
        print "As:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
            ExtMomBin[i], qData0[i]+qData1[i]/2.0, DataErr[i, 0]+DataErr[i, 1]/2.0)
        print "Aa:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
            ExtMomBin[i], qData1[i]/2.0, DataErr[i, 1]/2.0)
    else:
        print "Sum: {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
            ExtMomBin[i], qData0[i]+qData1[i], DataErr[i, 0]+DataErr[i, 1])


while True:

    time.sleep(SleepTime)

    for order in Order:
        for chan in Channel:

            files = os.listdir(folder)
            Num = 0
            Norm = 0
            Data0 = None
            DataList = []
            FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)

            for f in files:
                if re.match(FileName, f):
                    print "Loading ", f
                    Norm0 = -1
                    d = None
                    try:
                        with open(folder+f, "r") as file:
                            line0 = file.readline()
                            Step = int(line0.split(":")[-1])/1000000
                            # print "Step:", Step
                            line1 = file.readline()
                            # print line1
                            Norm0 = float(line1.split(":")[-1])
                            # print "Norm: ", Norm0
                            line3 = file.readline()
                            if AngleBin is None:
                                AngleBin = np.fromstring(
                                    line3.split(":")[1], sep=' ')
                                AngleBinSize = len(AngleBin)
                                print AngleBinSize
                            line4 = file.readline()
                            if ExtMomBin is None:
                                ExtMomBin = np.fromstring(
                                    line4.split(":")[1], sep=' ')
                                ExtMomBinSize = len(ExtMomBin)
                                ExtMomBin /= kF
                        # Num += 1
                        # print "Load data..."
                        d = np.loadtxt(folder+f)

                        if d is not None and Norm0 > 0:
                            if Data0 is None:
                                Data0 = d
                            else:
                                # Data0 = d
                                Data0 += d

                            Norm += Norm0

                            f = d.reshape(
                                (AngleBinSize, ExtMomBinSize, 2))/Norm0
                            DataList.append(AngleIntegation(f, 0))

                    # print "Norm", Norm

                    except:
                        print "fail to load ", folder+f

            if Norm > 0 and Data0 is not None:
                print "Total Weight: ", Data0[0]
                Data0 /= Norm
                Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize, 2))

                # print "Channel: ", chan
                if DataWithAngle.has_key((order, chan)):
                    DataWithAngle[(order, chan)] = DataWithAngle[(
                        order, chan)]*0.0+Data0*1.0
                else:
                    DataWithAngle[(order, chan)] = Data0

                Data[(order, chan)] = AngleIntegation(
                    DataWithAngle[(order, chan)], 0)

                # print np.array(DataList)
                DataErr[(order, chan)] = np.std(np.array(
                    DataList), axis=0)/np.sqrt(len(DataList))

                # print "err", np.std(np.array(DataList))

    if len(DataWithAngle) > 0:
        print "Write Weight file."
        for chan in Channel:
            with open("weight{0}.data".format(chan), "w") as file:
                for angle in range(AngleBinSize):
                    for qidx in range(ExtMomBinSize):
                        for Dir in range(2):
                            file.write("{0} ".format(
                                DataWithAngle[(0, chan)][angle, qidx, Dir]))

        with open("data.data", "a") as file:
            file.write("Dir: {0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(
                Data[(0, 1)][0, 0], Data[(0, 1)][0, 0], Data[(0, 2)][0, 0], Data[(0, 3)][0, 0]))
            file.write("Ex: {0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(
                Data[(0, 1)][0, 1], Data[(0, 1)][0, 1], Data[(0, 2)][0, 1], Data[(0, 3)][0, 1]))

        o = 0
        PrintInfo("T", Data[(o, 1)], DataErr[(o, 1)])
        PrintInfo("U", Data[(o, 2)], DataErr[(o, 2)])
        PrintInfo("S", Data[(o, 3)], DataErr[(o, 3)])
        qData = Data[(o, 1)]+Data[(o, 2)]+Data[(o, 3)]
        qDataErr = DataErr[(o, 1)]+DataErr[(o, 2)]+DataErr[(o, 3)]
        PrintInfo("Sum", qData, qDataErr)

    if Step >= TotalStep:
        print "End of Simulation!"
        sys.exit(0)
