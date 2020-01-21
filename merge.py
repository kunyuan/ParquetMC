import os
import sys
import re
import glob
import time
import numpy as np

SleepTime = 10

rs = None
Lambda = None
Beta = None
Charge2 = None
TotalStep = None
BetaStr = None
rsStr = None
ChargeStr = None
LambdaStr = None

with open("inlist", "r") as file:
    line = file.readline()
    para = line.split(" ")
    BetaStr = para[1]
    Beta = float(BetaStr)
    rsStr = para[2]
    rs = float(rsStr)
    LambdaStr = para[3]
    Lambda = float(LambdaStr)
    ChargeStr = para[4]
    Charge2 = float(ChargeStr)
    TotalStep = float(para[6])

print rs, Beta, Lambda, TotalStep

# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [0, ]

folder = "./Beta{0}_rs{1}_lambda{2}/".format(BetaStr, rsStr, LambdaStr)

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
                                Data0 += d

                            Norm += Norm0

                            f = d.reshape((AngleBinSize, ExtMomBinSize))/Norm0
                            DataList.append(AngleIntegation(f, 0))

                    # print "Norm", Norm

                    except:
                        print "fail to load ", folder+f

            if Norm > 0 and Data0 is not None:
                print "Total Weight: ", Data0[0]
                Data0 /= Norm
                Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize))

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
                        file.write("{0} ".format(
                            DataWithAngle[(0, chan)][angle, qidx]))

        with open("data.data", "a") as file:
            file.write("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(
                Data[(0, 1)][0], Data[(0, 1)][0], Data[(0, 2)][0], Data[(0, 3)][0]))

        qData = Data[(0, 1)]
        qData = 8.0*np.pi*Charge2/(ExtMomBin**2*kF**2+Lambda)-qData
        # print qData
        print "  Q/kF,    T,    Error"
        for i in range(len(qData)):
            print "{0:6.2f}, {1:10.6f}, {2:10.6f}".format(
                ExtMomBin[i], qData[i], DataErr[(0, 1)][i])

        qData = Data[(0, 2)]
        print "  Q/kF,    U,    Error"
        for i in range(len(qData)):
            print "{0:6.2f}, {1:10.6f}, {2:10.6f}".format(
                ExtMomBin[i], qData[i], DataErr[(0, 2)][i])

        qData = Data[(0, 3)]
        print "  Q/kF,    S,    Error"
        for i in range(len(qData)):
            print "{0:6.2f}, {1:10.6f}, {2:10.6f}".format(
                ExtMomBin[i], qData[i], DataErr[(0, 3)][i])

    if Step >= TotalStep:
        print "End of Simulation!"
        sys.exit(0)
