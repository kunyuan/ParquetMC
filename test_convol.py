import os
import sys
import re
import glob
import math
import numpy as np
from scipy import integrate as scp_int
from scipy import interpolate
#from pynufft import NUFFT_cpu, NUFFT_hsa
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from utility.IO import *
import grid




Order = [0, ]

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

def Extend_grid(grid):
    grid0=grid-grid[-1]-grid[0]
    grid1=grid[-1]+grid
    grid=np.concatenate([grid0,grid,grid1])
    return grid


def Extend_value(value):
    value0=-value
    value1=-value
    value=np.concatenate([value0,value,value1],axis=0)
    return value
    #print ("newgrid",grid)
    #print ("newfunction",value)

def Convol(ff,gg,Tau0,Size,swit):
    result=np.zeros((len(Tau0),Size))
    i=0
    for tau in Tau0:
        if(swit==-1):
            convol=scp_int.simps(ff(tau-Tau0)*gg(Tau0),Tau0)
        else:
            convol=scp_int.simps(ff(tau+Tau0)*gg(Tau0),Tau0)
        result[i]=convol
        i+=1
        #print (convol)
    return result
    

def CreateFolder(path):
    if not os.path.exists(path):
        os.system("mkdir "+path)

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


# with open("parameter","r") as file:
#     file.readline()
#     line0=file.readline()
#     print(line0.split(",")[-1])

#     #read Beta from F.dat
#     Beta = float(line0.split(",")[1])  
#     print ("Beta", Beta)
#     Temp=1/Beta

omega_c=10000000.0 #float(line0.split(",")[-1])  


Para = param()
Beta=Para.Beta
Temp=1/Beta
ExtMomBin=grid.FermiK()
ExtMomBin.build(1.0,3.0,8,1.0)

K=np.zeros(ExtMomBin.size)
for i in range (ExtMomBin.size):
    K[i]=ExtMomBin.grid[i]
#print (K)

#Order = range(0, Para.Order+1)
TauBin = grid.Tau()
TauBin.build(10.0, 16, 3.0)
#print(TauBin.str())
#print(TauBin.size)


Tau=np.zeros(TauBin.size)
for i in range (TauBin.size):
    Tau[i]=TauBin.grid[i]

N=256
Tau=np.arange(N)*Beta/N
print (Tau)
Tau0=Tau
FreqBin = (np.arange(TauBin.size)+0.5)*2*np.pi*Temp
FreqBinSize=len(FreqBin)

f_test=np.sin(np.pi*Tau/Beta)
f_test=np.tile(f_test,[3,1]).T
print (f_test.shape)
print ("oldfunction",f_test)

Tau=Extend_grid(Tau)
f_test=Extend_value(f_test)
f_test=f_test.T
#print ("newgrid",Tau)
print ("newfunction",f_test.shape)

ExtMomBinSize=3
ff=interpolate.interp1d(Tau,f_test)
#print ("interpolate",ff(Tau0))
gg=ff
result=Convol(ff,gg,Tau0,ExtMomBinSize,-1)
print(result)
f_test=Extend_value(result.T)
ff=interpolate.interp1d(Tau,f_test)
result=Convol(ff,gg,Tau0,ExtMomBinSize,1)

fig=plt.figure()
ax1=plt.axes()
ax1.plot(Tau0,result.T[0],'-')
ax1.plot(Tau0,np.sin(np.pi*Tau0)/4.0,'r-')
plt.show()

# #for order in Order:
#    # for chan in Channel:
# files = os.listdir(folder)
# Num = 0
# Norm = 0
# Data0 = None
# DataList = []
# #FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)
# FileName = "delta_pid0_verb.dat"
# FileName1="f.dat"




# #initialize delta
# gg=np.zeros((FreqBinSize,ExtMomBinSize))

# for i0 in range (FreqBinSize):
#     for j0 in range (ExtMomBinSize):
#         E=ExtMomBin[j0]*ExtMomBin[j0]-1.0
#         gg[i0][j0]=1.0 /(FreqBin[i0]*FreqBin[i0]+E*E)
    

# cut=-1
# delta=np.zeros(FreqBinSize*ExtMomBinSize)
# delta=delta.reshape((FreqBinSize,ExtMomBinSize))
# delta[0:cut,:]=1.0
# delta[cut:,:]=-0.1

# modulus=math.sqrt(np.tensordot(delta[0:cut,:],delta[0:cut,:],axes=([0,1],[0,1])))
# delta[0:cut,:]=delta[0:cut,:]/modulus


# IterationType=1
# loopcounter=0
# lamu=0
# shift=0.00
# lamu_sum=0.0
# modulus_dum=0.0

# #os set
# pid=0
# rootdir = os.getcwd()
# homedir = os.path.join(rootdir, "Data")
# CreateFolder(homedir)
# outfilepath = os.path.join(homedir, "outfile")
# CreateFolder(outfilepath)
# jobfilepath = os.path.join(homedir, "jobfile")
# CreateFolder(jobfilepath)
# seed = 1453
# outfile = os.path.join(outfilepath, "_out{0}".format(pid)) 
# jobfile = os.path.join(jobfilepath, "_job{0}.sh".format(pid))  # job files
# execute = "feyncalc.exe"




# f_freq=np.multiply(gg,delta)
# phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
# ff=phase_shift[:,np.newaxis]*np.fft.fft(f_freq,axis=0)/Beta
# F=2*ff.real


    
# order_num=3
# ll=2
# size0=len(d0[0])/order_num
# a=d0[1].reshape((order_num,size0))[ll].reshape((ExtMomBinSize,TauBinSize))
# d0=-d0[3].reshape((order_num,size0))[ll].reshape((ExtMomBinSize,TauBinSize))
# d0=np.transpose(d0)
# print d0
# ft_matrix=np.cos(np.tensordot(FreqBin,TauBin,axes=0))
# d=np.tensordot(ft_matrix,d0,axes=([1,0]))*Beta/TauBinSize

   
