#!/usr/bin/env python3
from utility.IO import *
import utility.fourier as fourier
import argparse
import time

parser = argparse.ArgumentParser("Specify some parameters.")
parser.add_argument("folder")
args = parser.parse_args()

folder = args.folder
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
    



if __name__ == "__main__":
    while True:
        Para = param(folder)
        Order = range(0, Para.Order+1)

        MaxFreq = 3000
        Freq = np.array(range(-MaxFreq, MaxFreq))
        phyFreq = (Freq*2.0+1.0)*np.pi/Para.Beta  # the physical frequency

        shape = (Para.Order+1, Para.MomGridSize, Para.TauGridSize)

        if not CheckFilesExist(folder, "sigma_pid[0-9]+.dat"):
            print("Files does not exist....")
            time.sleep(2)
            continue
        
        try:
            Data, Norm, Step, Grids = LoadFile(folder, "sigma_pid[0-9]+.dat", shape)
        except Exception as identifier:
            print("Error in loading files....")
            time.sleep(2)
            continue
        
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
        
        print("Mu=", Static[kFidx])
        Static -= Static[kFidx]  # subtract the self-energy shift

        print("Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr)))

        print("MomGrid idx at the Fermi surface:{0}, KF: {1}=={2}".format(kFidx,MomGrid[kFidx],Para.kF))
        # PlotSigmaW(Dynamic, MomGrid, kFidx, False)

        BareG = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)
        for i, q in enumerate(MomGrid):
            BareG[i, :] = 1.0/(1j*phyFreq+q*q-Static[i]-Para.EF)


        allorder, DynErr = Estimate(
            Data, Norm, lambda d: np.sum(d[1:Para.Order+1, ...], axis=0))
        SigmaW, _ = Fourier.SpectralT2W(allorder)

        s0, s1 = SigmaW[kFidx, MaxFreq-1], SigmaW[kFidx, MaxFreq]
        print("Z=", 1.0-(s1.imag-s0.imag)/(2.0*np.pi/Para.Beta))
        dMu = (s0.real+s1.real)/2.0
        print("Dynamic chemical shift: ", dMu)


        # dG_W = BareG
        dG_W = 1.0/(1.0/BareG+SigmaW-dMu)
        # dG_W = (SigmaW-dMu)*BareG*BareG/(1-(SigmaW-dMu)*BareG)

        dG_T, _ = Fourier.SpectralW2T(dG_W)
        dG_Tp = Fourier.naiveW2T(dG_W)
        # PlotSigmaT(dG_W, range(0, Para.MomGridSize, Para.MomGridSize/8), False)
        # print MomGrid[kFidx]/Para.kF
        # print BareG[kFidx, :]
        # PlotSigmaT(Static_MomTau, [kFidx], False)
        # PlotG(dG_T, kFidx, False)
        print("Maximum Error of \delta G: ", np.amax(abs(dG_T-dG_Tp)))

        with open(os.path.join(folder, "dispersion.data"), "w") as f:
            for k in range(Para.MomGridSize):
                f.write("{0} ".format(Static[k]))
            f.write("\n")

        with open(os.path.join(folder, "green.data"), "w") as f:
            for k in range(Para.MomGridSize):
                for t in range(Para.TauGridSize):
                    f.write("{0} ".format(dG_T[k, t].real))
            f.write("\n")

        with open(os.path.join(folder, "greenList.data"), "a+") as f:
            for k in range(Para.MomGridSize):
                for t in range(Para.TauGridSize):
                    f.write("{0} ".format(dG_T[k, t].real))
            f.write("\n")

        # with open("green_order.data", "w") as f:
        #     for k in range(Para.MomGridSize):
        #         for t in range(Para.TauGridSize):
        #             f.write("{0} ".format(dG_T[k, t].real))
        #     f.write("\n")

        time.sleep(2)

        print("\n")

        flag = np.array([step/1000000 >= Para.TotalStep for step in Step])
        print([step/1000000 for step in Step])
        print("Total Step: ", Para.TotalStep)
        if np.all(flag == True):
            print("End of Simulation!")
            sys.exit(0)
