import os
import sys
import re
import glob
import math
import numpy as np
#from pynufft import NUFFT_cpu, NUFFT_hsa
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



#rs = None
#Lambda = None
#Mass2 = None
#Beta = None
#Charge2 = None
#TotalStep = None
#BetaStr = None
#rsStr = None
#ChargeStr = None
#Mass2Str = None
#LambdaStr = None
#
#with open("inlist", "r") as file:
#    line = file.readline()
#    para = line.split(" ")
#    BetaStr = para[1]
#    Beta = float(BetaStr)
#    rsStr = para[2]
#    rs = float(rsStr)
#    Mass2Str = para[3]
#    Mass2 = float(Mass2Str)
#    LambdaStr = para[4]
#    Lambda = float(LambdaStr)
#    ChargeStr = para[5]
#    Charge2 = float(ChargeStr)
#    TotalStep = float(para[7])
#
#print rs, Beta, Mass2, Lambda, TotalStep
#Miu=1.0
# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [0, ]

#folder = "./Beta{0}_rs{1}_lambda{2}/".format(BetaStr, rsStr, Mass2Str)
folder = "./"

FreqBin=None
FreqBinSize=None
AngleBin = None
ExtMomBin = None
AngleBinSize = None
ExtMomBinSize = None
TauBin=None
TauBinSize=None
Data = {}  # key: (order, channel)
DataWithAngle = {}  # key: (order, channel)
DataErr = {}  # key: (order, channel)

Omega=3.0
g=2.0
kF = 1.0
Beta=1.0
Temp=1.0
#Nf = kF/2.0/np.pi**2
#Bubble = 0.0971916  # 3D, Beta=10, rs=1
#Step = None

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

#for order in Order:
   # for chan in Channel:
files = os.listdir(folder)
Num = 0
Norm = 0
Data0 = None
DataList = []
#FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)
FileName = "delta_pid0_verb.dat"
FileName1="delta0.dat"


for f in files:
    if re.match(FileName, f):
        #print ("Loading ")
        Norm0 = -1
        d = None
        #try: 
        with open(folder+f, "r") as file:
            #line1 = file.readline()
            # print (line1)
            #Norm0 = float(line1.split(":")[-1])
            # print ("Norm: ", Norm0)
            #file.open(folder+f, "r")
            d = np.transpose(np.loadtxt(folder+f,skiprows=2))

size0=len(d[0])/4
print size0

#phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
#ff=phase_shift[:,np.newaxis]*np.fft.fft(f_freq,axis=0)/Beta
#F=2*ff.real
fig=plt.figure()

ax=fig.add_subplot(111, projection='3d')
a=d[1].reshape((4,size0))[2]
print a
b=d[2].reshape((4,size0))[2]
print b
c=d[3].reshape((4,size0))[2]
print c
ax.scatter(a, b,c , c='k',marker='.')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()

