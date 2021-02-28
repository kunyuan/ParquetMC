#!/usr/bin/env python3
from utility.IO import *
import utility.fourier as fourier
import argparse
import time
# from utility.plot import *

parser = argparse.ArgumentParser("Specify some parameters.")
parser.add_argument("folder")
args = parser.parse_args()

folder = args.folder
print("------------------------ merge.py ------------------------------")
print("Folder to merge : " + folder)



def CheckConvergence(Para, Data):
    print(Para.Order)
    Order = range(0, Para.Order+1)
    print("    Q/kF,      Data,      Error")
    for o in Order[1:]:
        print("Order {0}".format(o))
        # sum all orders
        DataAllList = [np.sum(d[1:o+1, ...], axis=0) for d in Data]
        # average tau
        DataAllList = [np.average(d[:,:], axis=1) for d in DataAllList]
        
        Dat, Err = Estimate(DataAllList, Norm)
        print("tau average Sigma({0:2.1f}kF):, {1:10.6f}, {2:10.6f}".format(
            MomGrid[0], Dat[0], Err[0]))
    # sys.exit(0)

def CheckFilesExist(folder, fname):
    fexist = [bool(re.search(fname, f)) for f in getListOfFiles(folder)]
    return np.any(fexist)

def Sigma2VertexFolder(foldername):
    o = re.findall(r'(?<=_O)\d', foldername)[0]
    olast = str(int(o) + 1)
    fnew = foldername.replace("_O"+o, "_O*")
    fnew = re.sub(r'sigma_', 'F_', fnew)
    return fnew 



if __name__ == "__main__":
    Para = param(folder)
    selfFolder = os.path.join(folder, "selfconsistent")
    if not os.path.exists(selfFolder):
        os.system("mkdir " + selfFolder)
    
    Order = range(0, Para.Order+1)
    MaxFreq = 3000
    Freq = np.array(range(-MaxFreq, MaxFreq))
    phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency

    shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)

    if not CheckFilesExist(folder, "sigma_pid[0-9]+.dat"):
        print("Files does not exist....")
        time.sleep(2)
        sys.exit(0)
    
    try:
        Data, Norm, Step, Grids = LoadFile(folder, "sigma_pid[0-9]+.dat", shape)
    except Exception as identifier:
        print("Error in loading files....")
        time.sleep(2)
        sys.exit(0)
    
    TauGrid = Grids["TauGrid"]
    MomGrid = Grids["KGrid"]

    # Data.shape : (pid, order+1, MomGrid, TauGrid)
    
    CheckConvergence(Para, Data)

    Fourier = fourier.fourier(TauGrid, phyFreq, Para.Beta)
    Fourier.InitializeKernel(100.0, 1024, "Fermi", 1.0e-13)

    # first order is a constant of tau
    # Static(MomGrid) : order 1
    Static, StaticErr = Estimate(
        Data, Norm, lambda d: np.average(d[1, :, :], axis=1))
    # Dynamic(MomGrid, TauGrid) : order 2-n
    Dynamic, DynErr = Estimate(
        Data, Norm, lambda d: np.sum(d[2:Para.Order+1, ...], axis=0))
    
    arr = np.amin(abs(MomGrid-Para.kF))
    kFidx = np.where(abs(arr - abs(MomGrid-Para.kF)) < 1.0e-20)[0][0]
    
    print("Mu=", Static[kFidx], " +- ", StaticErr[kFidx])
    Static -= Static[kFidx]  # subtract the self-energy shift

    # def PlotStatic():
    #     fig, ax1 = plt.subplots()
        
    #     ax1.plot(MomGrid/Para.kF, 2*Static, "r-", label="Static Fock")
    #     ax1.set_ylabel('$\Sigma_{Fock}$')
    #     ax1.set_xlabel("$k/k_F$")

    #     plt.legend(loc=4, frameon=False, fontsize=size)
    #     plt.show()
    # PlotStatic()
    # sys.exit(0)

    print("Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr)))

    print("MomGrid idx at the Fermi surface:{0}, KF: {1}=={2}".format(kFidx,MomGrid[kFidx],Para.kF))
    # PlotSigmaW(Dynamic, MomGrid, kFidx, False)
    
    fname = "para.data"
    with open(os.path.join(selfFolder,fname), "w") as f:
        f.write("{0}  {1}  {2}\n".format(Para.TauGridSize, Para.MomGridSize, Para.MaxExtMomKF))
        f.write("#TauGrid, MomGrid, MaxExtMom(*kF)")


    fname = "dispersion_order{0}.data".format(Para.Order)
    with open(os.path.join(selfFolder,fname), "w") as f:
        for k in range(Para.MomGridSize):
            f.write("{0} ".format(Static[k]))
        f.write("\n")

    for o in range(Para.Order, Para.Order+1):
        Dynamic, DynErr = Estimate(Data, Norm, lambda d: np.sum(d[2:o+1, ...], axis=0))
        SigmaW, _ = Fourier.SpectralT2W(Dynamic)
        SigmaWErr, _ = Fourier.SpectralT2W(DynErr)

        # with open("Dynamic.data", "w") as f:
        #     for i1 in range(len(Dynamic)):
        #         for j1 in range(len(Dynamic[i1])):
        #             f.write("{0}  ".format(Dynamic[i1, j1]))
        #         f.write("\n")
        # with open("Static.data", "w") as f:
        #     for v in Static:
        #         f.write("{0}  ".format(v))
        # with open("MomGrid.data", "w") as f:
        #     for v in MomGrid:
        #         f.write("{0}  ".format(v))
        # with open("TauGrid.data", "w") as f:
        #     for i1 in TauGrid:
        #         f.write("{0}  ".format(i1))

        
        # ================ Z, dMu and m* ================
        s0, s1 = SigmaW[kFidx, MaxFreq-1], SigmaW[kFidx, MaxFreq]
        Z = 1.0-(s1.imag-s0.imag)/(2.0*np.pi/Para.Beta)
        print("order={0}  Z={1}".format(o, Z) )

        dMu = (s0.real+s1.real)/2.0
        dMuErr = (SigmaWErr[kFidx, MaxFreq-1] + SigmaWErr[kFidx, MaxFreq])/2.0
        print("Dynamic chemical shift: ", dMu, "+-", dMuErr)

        sigmaW0 = [0.0] * len(MomGrid)
        for i1, k1 in enumerate(MomGrid):
            dynamicW0 = sum([(Dynamic[i1,j]+Dynamic[i1,j+1])/2.0 * (TauGrid[j+1]-TauGrid[j]) for j in range(len(TauGrid)-1)])
            sigmaW0[i1] = Static[i1] + dynamicW0

        rangeL, rangeR = kFidx - 3, kFidx + 3 + 1
        coeff = np.polyfit(MomGrid[rangeL:rangeR], sigmaW0[rangeL:rangeR], 1)

        m = 1.0/2.0
        mStar = m/(Z*(1 + coeff[0]*m/Para.kF))
        print("m*={0}, m*/m={1}".format(mStar, mStar/m))

        fname = "ZandM_order{0}.data".format(o)
        with open(os.path.join(selfFolder, fname), "w") as f:
            f.write("{0}    {1}".format(Z, mStar))

        # ================ Z, dMu and m* ================

        BareGW = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)
        for i, q in enumerate(MomGrid):
            BareGW[i, :] = 1.0/(1j*phyFreq + ( (q*q-Para.EF) + Static[i] ) )


        BoldGW = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)
        for i, q in enumerate(MomGrid):
            for j, w in enumerate(phyFreq):
                BoldGW[i, j] = 1.0/(1j*w + ( (q*q-Para.EF) + Static[i] + (SigmaW[i,j]-dMu) ) )

        dG_W = BoldGW - BareGW

        dG_T, _ = Fourier.SpectralW2T(dG_W)

        fname = "green_order{0}.data".format(o)
        with open(os.path.join(selfFolder, fname), "w") as f:
            for k in range(Para.MomGridSize):
                for t in range(Para.TauGridSize):
                    f.write("{0} ".format(dG_T[k, t].real))
            f.write("\n")


    cpFolder = Sigma2VertexFolder(selfFolder)
    os.system("cp {0}/*  {1}/".format(selfFolder, cpFolder))
    print("--------------------- End merge.py ---------------------------")
    sys.exit(0)
