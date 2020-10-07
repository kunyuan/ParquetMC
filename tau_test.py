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
folder = "./Data/"

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
#Nf = kF/2.0/np.pi**2
#Bubble = 0.0971916  # 3D, Beta=10, rs=1
#Step = None

def Extend_grid(grid):
    grid0=grid-grid[-1]-grid[0]
    grid1=grid[-1]+grid
    grid=np.concatenate([grid0,grid,grid1])
    return grid


def Extend_value(value):
    value0=-value
    value1=-value
    value=np.concatenate([value0,value,value1],axis=1)
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
    return result.T
    


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


# with open("Data/grid.data","r") as file:
#     file.readline()
#     line2 = file.readline()
#     if TauBin is None:
#         TauBin = np.fromstring(
#             line2, sep=' ')
#         TauBinSize = len(TauBin)
#         print (TauBinSize)
#     #line3 = file.readline()
#     #if AngleBin is None:
#     #    AngleBin = np.fromstring(
#       line3.split(":")[1], sep=' ')
#     #    AngleBinSize = len(AngleBin)
#     #    print (AngleBinSize)
#     file.readline()
#     file.readline()
#     line4 = file.readline()
#     if ExtMomBin is None:
#         ExtMomBin = np.fromstring(
#             line4, sep=' ')
#         ExtMomBinSize = len(ExtMomBin)
#         ExtMomBin /= kF   
# print ExtMomBin

Para = param()
Beta=Para.Beta
Temp=1/Beta
K=grid.FermiK()
K.build(1.0,3.0,8,1.0) #kf,maxk,size,scale
Ta=grid.Tau()
Ta.build(10.0, 16, 3.0) #Beta,size,scale

#FreqBin = (np.arange(len(TauBin))+0.5)*2*np.pi*Temp
#FreqBinSize=len(FreqBin)
ExtMomBin=np.zeros(K.size)
for i in range (K.size):
    ExtMomBin[i]=K.grid[i]

TauBin=np.zeros(Ta.size)
for i in range (Ta.size):
    TauBin[i]=Ta.grid[i]

N=256
TauBin=np.arange(N)*Beta/N
TauBinSize=len(TauBin)
ExtMomBinSize=len(ExtMomBin)

TauExtend=Extend_grid(TauBin)

omega_c=10000000.0 #float(line0.split(",")[-1])  

#for order in Order:
   # for chan in Channel:
files = os.listdir(folder)
Num = 0
Norm = 0
Data0 = None
DataList = []
#FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)
FileName = "delta_pid0_verb.dat"
FileName1="f.dat"



#initialize F


cut=-1#np.searchsorted(TauBin,100000,side='right') 

# if(os.path.exists(folder+FileName1) and os.path.getsize(folder+FileName1)) > 0:
#     F = np.loadtxt(folder+FileName1)
#     F=F.reshape((TauBinSize,ExtMomBinSize))


# else:
#     F=np.zeros(TauBinSize*ExtMomBinSize)
#     F=F.reshape((TauBinSize,ExtMomBinSize))
#     F[:,:]=1.0


#phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
#dd=phase_shift[:,np.newaxis]*np.fft.fft(F,axis=0)/Beta
#ff=np.multiply(gg_tau,dd)
#F=2*ff.real



#modulus=math.sqrt(np.tensordot(F[0:cut,:],F[0:cut,:],axes=([0,1],[0,1])))
#F[0:cut,:]=F[0:cut,:]/modulus


IterationType=1
loopcounter=0
lamu=0
shift=0.00
lamu_sum=0.0
modulus_dum=0.0

#os set
pid=0
rootdir = os.getcwd()
homedir = os.path.join(rootdir, "Data")
CreateFolder(homedir)
outfilepath = os.path.join(homedir, "outfile")
CreateFolder(outfilepath)
seed = 1453
outfile = os.path.join(outfilepath, "_out{0}".format(pid)) 
execute = "feyncalc.exe"


for loopcounter in range(1):


    #phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
    #dd=phase_shift[:,np.newaxis]*np.fft.fft(F,axis=0)/Beta
    #ff=np.multiply(gg_tau,dd)
    #F=2*ff.real


    # with open("./Data/f.dat","w") as file:
    #     for i in range(TauBinSize):
    #         file.write("{0} ".format(TauBin[i]))
    #     file.write("\n")
    #     for i in range(TauBinSize):
    #         for k in range(ExtMomBinSize):
    #             file.write("{0}\t".format(F[i][k]))
    
    #if(loopcounter%5==0):
     #  IterationType=(IterationType+1)%2
    #print IterationType
    #os.chdir(homedir)
    #os.system("./{0} {1} {2} > {3} &".format(execute, pid, seed, outfile))
    #os.chdir("..")

     #myCmd='python phonon.py>out.txt'

    # for f in files:
    #     if re.match(FileName, f):
    #         #print ("Loading ")
    #         Norm0 = -1
    #         d0 = None
    #         #try: 
    #         with open(folder+f, "r") as file:
    #             #line1 = file.readline()
    #             # print (line1)
    #             #Norm0 = float(line1.split(":")[-1])
    #             # print ("Norm: ", Norm0)
    #             #file.open(folder+f, "r")
    #             d0 = np.transpose(np.loadtxt(folder+f))
    #             #fig=plt.figure()
    #             #plt.axvline(omega_c)
    #             #plt.plot(FreqBin,F[:,0],'ko')
    #             #plt.plot(ExtMomBin,F[0,:],'r.')
    #             #ax=fig.add_subplot(111, projection='3d')
    #             #a=FreqBin[:,np.newaxis]+ExtMomBin*0
    #             #a=a.reshape(TauBinSize*ExtMomBinSize)
    #             #print a
    #             #b=np.transpose(ExtMomBin[:,np.newaxis]+FreqBin*0)
    #             #b=b.reshape(TauBinSize*ExtMomBinSize)
    #             #print b
    #             #c=d.reshape(TauBinSize*ExtMomBinSize)
    #             #modulus=math.sqrt(np.tensordot(c,c,axes=(0,0)))
    #             #c=c/modulus
    #             #ax.scatter(a, b,c , c='k',marker='.')
    #             #ax.set_xlabel('X Label')
    #             #ax.set_ylabel('Y Label')
    #             #ax.set_zlabel('Z Label')

    #             #plt.savefig('d_q.png')
    #             #plt.close()
    #             #plt.figure()
    #             #plt.axvline(omega_c)
    #             #plt.plot(FreqBin,F[:,0],'ko')
    #             #plt.plot(ExtMomBin,F[0,:],'r.')
    #             #plt.savefig('F_new.pdf')
    #             #plt.close()
    #             #if d is not None and Norm0 > 0:
    #                 #if Data0 is None:
    #                     #Data0 = d
    #                 #else:
    #                     # Data0 = d
    #                 #   Data0 += 
            
   
#       print (TauBin)
    # order_num=3
    # ll=2
    # size0=len(d0[0])/order_num
    # d=-d0[3].reshape((order_num,size0))[ll].reshape((TauBinSize,ExtMomBinSize))
    # print d
    d=np.sin(np.pi*TauBin/Beta)
    d=np.tile(d,[ExtMomBinSize,1])
    print (d.shape)
    #do twice convolution to convert delta to F
    d=Extend_value(d)
    d_int=interpolate.interp1d(TauExtend,d)
    g=d_int
    middle=Convol(d_int,g,TauBin,ExtMomBinSize,-1)
    middle=Extend_value(middle)
    d_int=interpolate.interp1d(TauExtend,middle)
    middle=Convol(d_int,g,TauBin,ExtMomBinSize,1)

    fig=plt.figure()
    ax1=plt.axes()
    ax1.plot(TauBin,middle[0],'-')
    ax1.plot(TauBin,np.sin(np.pi*TauBin)/4.0,'r-')
    plt.show()
    # F0=0.0#-g/4.0/np.pi/np.pi*sum(ExtMomBin*F[0])/ExtMomBinSize*ExtMomBin[-1]
    # taudep= np.cosh(Omega*(0.5*Beta-np.abs(TauBin))) * np.tensordot(ExtMomBin,F,axes=([0,1]))/ExtMomBinSize*ExtMomBin[-1]
    # taudep=taudep*g*Omega/8.0/np.pi/np.pi/np.sinh(0.5*Beta*Omega)+F0
    # d+=np.tensordot(taudep,ExtMomBin,axes=0)
    
    #test=np.sin((Beta*FreqBin-np.pi)/TauBinSize)
# print (test)
# dd[:,0,0,0]=dd[:,0,0,0]+test
    #separate F in to high and low frequency 

    
    # if(IterationType==0):
    #     F[cut:,:]=d[cut:,:]
    # elif(IterationType==1): 
    #     d[0:cut,:]=d[0:cut,:]+(shift+0.9*modulus_dum)*F[0:cut,:]
    #     lamu=np.tensordot(d[0:cut,:],F[0:cut,:],axes=([0,1],[0,1]))
    #     lamu_sum=lamu_sum*0.9+lamu-shift-modulus_dum*0.9
    #     print lamu-shift-0.9*modulus_dum
    #     print lamu_sum/10.0
    #     modulus_dum=math.sqrt(np.tensordot(d[0:cut,:],d[0:cut,:],axes=([0,1],[0,1])))
    #     F[0:cut,:]=d[0:cut,:]/modulus_dum  

    # with open("lamu.txt","a+") as file:
    # 	file.write("{0} \n".format(lamu_sum/10.0))




    #fig=plt.figure()
    #plt.axvline(omega_c)
    #plt.plot(FreqBin,F[:,0],'ko')
    #plt.plot(ExtMomBin,F[0,:],'r.')
    #ax=fig.add_subplot(111, projection='3d')
    #a=FreqBin[:,np.newaxis]+ExtMomBin*0
    #a=a.reshape(TauBinSize*ExtMomBinSize)
    #print a
    #b=np.transpose(ExtMomBin[:,np.newaxis]+FreqBin*0)
    #b=b.reshape(TauBinSize*ExtMomBinSize)
    #print b
    #c=F.reshape(TauBinSize*ExtMomBinSize)
    #ax.scatter(a, b,c , c='k',marker='.')
    #ax.set_xlabel('X Label')
    #ax.set_ylabel('Y Label')
    #ax.set_zlabel('Z Label')
    
    #F.reshape((TauBinSize,ExtMomBinSize))

    #plt.savefig('F_q.png')
    #plt.close()
    
    #cut1_num=10.0
    #cut2_num=3.0
    #cut1=np.searchsorted(FreqBin,cut1_num,side='right')
    #cut2=np.searchsorted(ExtMomBin,cut2_num,side='right')
    #for i in range(cut1):
    #    for j in range(cut2):
    #        if(j==cut2-1):
    #            step2=cut2_num-ExtMomBin[j]
    #        else:
    #            step2=ExtMomBin[j+1]-ExtMomBin[j]
    #        step1=2*np.pi/Beta
    #        sum0+=F[i][j]*step2*step1
    #plt.show()



  #  fig=plt.figure()
  #  plt.plot(FreqBin,F[:,-1],'k.')
  #  plt.savefig('F_q0.png')
  #  plt.close()

  #  fig=plt.figure()
  #  plt.plot(ExtMomBin,F[0,:],'k.')
  #  plt.savefig('F_q1.png')
  #  plt.close()

    ##Fourier Transformation of F*GG 

# om=np.random.randn(1512,1) #None-Uniform-Fourier-Transform
# Nd=(256,)
# Kd=(512,)
# Jd=(6,)
    
# time_data = np.zeros(256, )
# time_data[96:128+32] = 1.0

# NufftObj=NUFFT_cpu()
# NufftObj.plan(om,Nd,Kd,Jd)
# ff=NufftObj.forward(time_data)
        
# Fourier=np.exp(1j*2*np.pi*np.tensordot(FreqBin,TauBin,0))
# print Fourier.shape
# print Fourier
#  ff=np.tensordot(Fourier,dd,axes=([0,0]))
                #  Norm += Norm0
                #   dd = d.reshape(
                #       (AngleBinSize, ExtMomBinSize, 2))/Norm0
                #   dd = AngleIntegation(dd, 0)
                #   DataList.append(SpinMapping(dd))
                        
        # except:
            #    print "fail to load ", folder+f
            
    #if Norm > 0 and Data0 is not None:
    #        print "Total Weight: ", Data0[0]
    #        Data0 /= Norm
    #        Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize, 2)

   
