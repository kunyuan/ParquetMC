import os
import sys
import re
import glob
import math
import numpy as np
#from pynufft import NUFFT_cpu, NUFFT_hsa
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from utility.IO import *
import grid
import psutil


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
#Nf = kF/2.0/np.pi**2
#Bubble = 0.0971916  # 3D, Beta=10, rs=1
#Step = None
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


with open("parameter","r") as file:
    file.readline()
    line0=file.readline()
    print line0.split(",")[-1]

    #read Beta from F.dat
    Beta = float(line0.split(",")[1])  
    print ("Beta", Beta)
    Temp=1/Beta

omega_c=10000000.0 #float(line0.split(",")[-1])  


Para = param()
ExtMomBin=BuildMomGrid(Para.MaxExtMom, MomGridSize)

#Order = range(0, Para.Order+1)
TauBin = Grid.TauGrid
TauBinSize=len(TauBin)
FreqBin = (np.arange(TauBinSize)+0.5)*2*np.pi*Temp
FreqBinSize=len(FreqBin)

ExtMomBinSize=len(ExtMomBin)



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




#initialize delta
gg=np.zeros((FreqBinSize,ExtMomBinSize))

for i0 in range (FreqBinSize):
    for j0 in range (ExtMomBinSize):
        E=ExtMomBin[j0]*ExtMomBin[j0]-1.0
        gg[i0][j0]=1.0 /(FreqBin[i0]*FreqBin[i0]+E*E)
    

cut=-1#np.searchsorted(FreqBin,omega_c,side='right') 

#if(os.path.exists(folder+FileName1) and os.path.getsize(folder+FileName1)) > 0:
#    delta = np.loadtxt(folder+FileName1)
#    delta=delta.reshape((FreqBinSize,ExtMomBinSize))


#else:
delta=np.zeros(FreqBinSize*ExtMomBinSize)
delta=delta.reshape((FreqBinSize,ExtMomBinSize))
delta[0:cut,:]=1.0
delta[cut:,:]=-0.1
# print delta
# f_freq=np.multiply(gg,delta)
# phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
# ff=phase_shift[:,np.newaxis]*np.fft.fft(f_freq,axis=0)/Beta
# F=2*ff.real

# #phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
# #dd=phase_shift[:,np.newaxis]*np.fft.fft(delta,axis=0)/Beta
# #ff=np.multiply(gg_tau,dd)
# #F=2*ff.real


# delta0=0.0#-g/4.0/np.pi/np.pi*sum(ExtMomBin*F[0])/ExtMomBinSize*ExtMomBin[-1]
# taudep= np.cosh(Omega*(0.5*Beta-np.abs(TauBin))) * np.tensordot(ExtMomBin,F,axes=([0,1]))/ExtMomBinSize*ExtMomBin[-1]
# delta1=np.tensordot(taudep,np.cos(np.tensordot(FreqBin,TauBin,axes=0)),axes=([0,1]))/TauBinSize*Beta
# delta1=delta1*g*Omega/8.0/np.pi/np.pi/np.sinh(0.5*Beta*Omega)+delta0
# delta+=np.tensordot(delta1,ExtMomBin,axes=0)

modulus=math.sqrt(np.tensordot(delta[0:cut,:],delta[0:cut,:],axes=([0,1],[0,1])))
delta[0:cut,:]=delta[0:cut,:]/modulus


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
jobfilepath = os.path.join(homedir, "jobfile")
CreateFolder(jobfilepath)
seed = 1453
outfile = os.path.join(outfilepath, "_out{0}".format(pid)) 
jobfile = os.path.join(jobfilepath, "_job{0}.sh".format(pid))  # job files
execute = "feyncalc.exe"
os.system("cp {0} {1}".format(execute, homedir))
os.system("cp {0} {1}".format("parameter", homedir))

os.system("python grid.py")
os.system("cp {0} {1}".format("grid.data", homedir))

os.chdir(homedir)


for loopcounter in range(1):

    f_freq=np.multiply(gg,delta)
    phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
    ff=phase_shift[:,np.newaxis]*np.fft.fft(f_freq,axis=0)/Beta
    F=2*ff.real


    with open("./f.dat","w") as file:
        for i in range(TauBinSize):
            file.write("{0} ".format(TauBin[i]))
        file.write("\n")
        for i in range(TauBinSize):
            for k in range(ExtMomBinSize):
                if(loopcounter==0):
                    file.write("{0}\t".format( np.exp(-ExtMomBin[k]**2)*(Para.Beta-2*TauBin[i]) ))
                else:
                    file.write("{0}\t".format(F[i][k]))


    #if(loopcounter%5==0):
     #  IterationType=(IterationType+1)%2
    print IterationType
    
    

    
    os.system("./{0} {1} {2} >{3}".format(execute, pid, seed, outfile))
    
    # for proc in psutil.process_iter():
    #     try:
    #         # this returns the list of opened files by the current process
    #         flist = proc.open_files()
            
    # # This catches a race condition where a process ends
    # # before we can examine its files    
    #     except psutil.NoSuchProcess as err:
    #         print("****",err) 



    for f in files:
        if re.match(FileName, f):
            #print ("Loading ")
            Norm0 = -1
            d0 = None
            #try: 
            with open(folder+f, "r") as file:
                #line1 = file.readline()
                # print (line1)
                #Norm0 = float(line1.split(":")[-1])
                # print ("Norm: ", Norm0)
                #file.open(folder+f, "r")
                d0 = np.transpose(np.loadtxt(folder+f))
                #fig=plt.figure()
                #plt.axvline(omega_c)
                #plt.plot(FreqBin,delta[:,0],'ko')
                #plt.plot(ExtMomBin,delta[0,:],'r.')
                #ax=fig.add_subplot(111, projection='3d')
                #a=FreqBin[:,np.newaxis]+ExtMomBin*0
                #a=a.reshape(FreqBinSize*ExtMomBinSize)
                #print a
                #b=np.transpose(ExtMomBin[:,np.newaxis]+FreqBin*0)
                #b=b.reshape(FreqBinSize*ExtMomBinSize)
                #print b
                #c=d.reshape(FreqBinSize*ExtMomBinSize)
                #modulus=math.sqrt(np.tensordot(c,c,axes=(0,0)))
                #c=c/modulus
                #ax.scatter(a, b,c , c='k',marker='.')
                #ax.set_xlabel('X Label')
                #ax.set_ylabel('Y Label')
                #ax.set_zlabel('Z Label')

                #plt.savefig('d_q.png')
                #plt.close()
                #plt.figure()python
                #plt.axvline(omega_c)
                #plt.plot(FreqBin,delta[:,0],'ko')
                #plt.plot(ExtMomBin,delta[0,:],'r.')
                #plt.savefig('delta_new.pdf')
                #plt.close()
                #if d is not None and Norm0 > 0:
                    #if Data0 is None:
                        #Data0 = d
                    #else:
                        # Data0 = d
                    #   Data0 += 
            
   
#       print (TauBin)
    order_num=3
    ll=2
    size0=len(d0[0])/order_num
    a=d0[1].reshape((order_num,size0))[ll].reshape((ExtMomBinSize,TauBinSize))
    d0=-d0[3].reshape((order_num,size0))[ll].reshape((ExtMomBinSize,TauBinSize))
    d0=np.transpose(d0)
    print d0
    ft_matrix=np.cos(np.tensordot(FreqBin,TauBin,axes=0))
    d=np.tensordot(ft_matrix,d0,axes=([1,0]))*Beta/TauBinSize
    #d=d0
    # FileName1="/home/wangtao/Parquet_test3/w0.dat"
    # d2=np.transpose(np.loadtxt(FileName1))
    
    # fig=plt.figure()
    # ax1=plt.axes()
    # aa=a.T[0]
    # dd=d[0]
    # #print aa
    # print aa[32],dd[32]
    # ax1.plot(aa,-dd,label="new")

    # ax1.plot(d2[1],-d2[2],label="old")

    # ax1.legend(loc=[0.7,0.7], shadow=False,fontsize=10)


    # plt.show()


    # delta0=0.0#-g/4.0/np.pi/np.pi*sum(ExtMomBin*F[0])/ExtMomBinSize*ExtMomBin[-1]
    # taudep= np.cosh(Omega*(0.5*Beta-np.abs(TauBin))) * np.tensordot(ExtMomBin,F,axes=([0,1]))/ExtMomBinSize*ExtMomBin[-1]
    # delta1=np.tensordot(taudep,np.cos(np.tensordot(FreqBin,TauBin,axes=0)),axes=([0,1]))/TauBinSize*Beta
    # delta1=delta1*g*Omega/8.0/np.pi/np.pi/np.sinh(0.5*Beta*Omega)+delta0
    # d+=np.tensordot(delta1,ExtMomBin,axes=0)
    
    #test=np.sin((Beta*FreqBin-np.pi)/FreqBinSize)
# print (test)
# dd[:,0,0,0]=dd[:,0,0,0]+test
    #separate delta in to high and low frequency 


    lamu=np.tensordot(d[0:cut,:],delta[0:cut,:],axes=([0,1],[0,1]))
    print lamu
    modulus_dum=math.sqrt(np.tensordot(d[0:cut,:],d[0:cut,:],axes=([0,1],[0,1])))
    print ("modulus:",modulus_dum)
    delta[0:cut,:]=d[0:cut,:]/modulus_dum
    
    #print (FreqBin)
    # if(IterationType==0):
    #     delta[cut:,:]=d[cut:,:]
    # elif(IterationType==1): 
    #     d[0:cut,:]=d[0:cut,:]+(shift+0.9*modulus_dum)*delta[0:cut,:]
    #     lamu=np.tensordot(d[0:cut,:],delta[0:cut,:],axes=([0,1],[0,1]))
    #     lamu_sum=lamu_sum*0.9+lamu-shift-modulus_dum*0.9
    #     print lamu-shift-0.9*modulus_dum
    #     print lamu_sum/10.0
    #     modulus_dum=math.sqrt(np.tensordot(d[0:cut,:],d[0:cut,:],axes=([0,1],[0,1])))
    #     print ("modulus:",modulus)
    #     delta[0:cut,:]=d[0:cut,:]/modulus_dum  
    
    # with open("lamu.txt","a+") as file:
    # 	file.write("{0} \n".format(lamu_sum/10.0))




    #fig=plt.figure()
    #plt.axvline(omega_c)
    #plt.plot(FreqBin,delta[:,0],'ko')
    #plt.plot(ExtMomBin,delta[0,:],'r.')
    #ax=fig.add_subplot(111, projection='3d')
    #a=FreqBin[:,np.newaxis]+ExtMomBin*0
    #a=a.reshape(FreqBinSize*ExtMomBinSize)
    #print a
    #b=np.transpose(ExtMomBin[:,np.newaxis]+FreqBin*0)
    #b=b.reshape(FreqBinSize*ExtMomBinSize)
    #print b
    #c=delta.reshape(FreqBinSize*ExtMomBinSize)
    #ax.scatter(a, b,c , c='k',marker='.')
    #ax.set_xlabel('X Label')
    #ax.set_ylabel('Y Label')
    #ax.set_zlabel('Z Label')
    
    #
    #delta.reshape((FreqBinSize,ExtMomBinSize))

    #plt.savefig('delta_q.png')
    #plt.close()
    
    #cut1_num=10.0
    #cut2_num=3.0
    #cut1=np.searchsorted(FreqBin,cut1_num,side='right')
    #cut2=np.searchsorted(ExtMomBin,cut2_num,side='right')
    # cutf=np.searchsorted(ExtMomBin,1.0,side='right')
    # sum0=1.0
    # #for i in range(cut1):
    # #    for j in range(cut2):
    # #        if(j==cut2-1):
    # #            step2=cut2_num-ExtMomBin[j]
    # #        else:
    # #            step2=ExtMomBin[j+1]-ExtMomBin[j]
    # #        step1=2*np.pi/Beta
    # #        sum0+=delta[i][j]*step2*step1
    # #plt.show()

    # with open("qmax.dat","w") as file:
    #     for i in range(FreqBinSize):
    #         file.write("{0} ".format(FreqBin[i]))
    #     file.write("\n")
    #     for i in range(FreqBinSize):
    #         file.write("{0} ".format(delta[i][-1]/sum0))
    
    # with open("q0.dat","w") as file:
    #     for i in range(FreqBinSize):
    #         file.write("{0} ".format(FreqBin[i]))
    #     file.write("\n")
    #     for i in range(FreqBinSize):
    #         file.write("{0} ".format(delta[i][0]/sum0))

    # with open("w0.dat","w") as file:
    #     for i in range(ExtMomBinSize):
    #         file.write("{0} ".format(ExtMomBin[i]))
    #     file.write("\n")
    #     for i in range(ExtMomBinSize):
    #        file.write("{0} ".format(delta[0][i]/sum0))

    # with open("Fermi.dat","w") as file:
    #     for i in range(FreqBinSize):
    #         file.write("{0} ".format(FreqBin[i]))
    #     file.write("\n")
    #     for i in range(FreqBinSize):
    #         file.write("{0} ".format(delta[i][cutf]/sum0))


  #  fig=plt.figure()
  #  plt.plot(FreqBin,delta[:,-1],'k.')
  #  plt.savefig('delta_q0.png')
  #  plt.close()

  #  fig=plt.figure()
  #  plt.plot(ExtMomBin,delta[0,:],'k.')
  #  plt.savefig('delta_q1.png')
  #  plt.close()

    ##Fourier Transformation of delta*GG 

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

   
