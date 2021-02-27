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


def AnaliticFock(k, l):
    kF = Para.kF
    # Analytic Fock energy
    y = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((k-kF)/l)-l/kF*np.arctan((k+kF)/l) -
                      (l*l-k*k+kF*kF)/4.0/k/kF*np.log((l*l+(k-kF)**2)/(l*l+(k+kF)**2)))

    k = kF
    Mu = 2.0*kF/np.pi*(1.0+l/kF*np.arctan((k-kF)/l)-l/kF*np.arctan((k+kF)/l) -
                       (l*l-k*k+kF*kF)/4.0/k/kF*np.log((l*l+(k-kF)**2)/(l*l+(k+kF)**2)))
    return y-Mu


def PlotStatic():
    fig, ax1 = plt.subplots()
    fock = [AnaliticFock(k, Para.Lambda) for k in MomGrid]

    print(fock)
    ax1.plot(MomGrid/Para.kF, fock, "r-", label="Analitical Fock")
    ax1.plot(MomGrid/Para.kF, 1.5*Static, "b-", label="Static")

    plt.legend(loc=1, frameon=False, fontsize=size)
    plt.show()


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

        # print(Static)
        # PlotStatic()

        print("Maximum Error of Dynamic Sigma: ", np.amax(abs(DynErr)))

        print("MomGrid idx at the Fermi surface:{0}, KF: {1}=={2}".format(kFidx,MomGrid[kFidx],Para.kF))
        # PlotSigmaW(Dynamic, MomGrid, kFidx, False)

        with open(os.path.join(selfFolder,"dispersion.data"), "w") as f:
            for k in range(Para.MomGridSize):
                f.write("{0} ".format(Static[k]))
            f.write("\n")

        for o in range(2, Para.Order+1):
            Dynamic, DynErr = Estimate(Data, Norm, lambda d: np.sum(d[2:o+1, ...], axis=0))
            SigmaW, _ = Fourier.SpectralT2W(Dynamic)
            
            s0, s1 = SigmaW[kFidx, MaxFreq-1], SigmaW[kFidx, MaxFreq]
            Z = 1.0-(s1.imag-s0.imag)/(2.0*np.pi/Para.Beta)
            print("order={0}\n Z={1}".format(o, Z) )
            dMu = (s0.real+s1.real)/2.0
            print("Dynamic chemical shift: ", dMu)

            BareGW = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)
            for i, q in enumerate(MomGrid):
                # BareGW[i, :] = Z/(1j*phyFreq + (q*q-Para.EF) +Static[i] )
                BareGW[i, :] = 1.0/(1j*phyFreq + (q*q-Para.EF) + Static[i] )


            BoldGW = np.zeros((Para.MomGridSize, len(phyFreq)), dtype=complex)
            for i, q in enumerate(MomGrid):
                for j, w in enumerate(phyFreq):
                    # BoldGW[i, j] = Z/( 1j*w + (q*q-Para.EF) + Static[i] + (SigmaW[i,j]-dMu) )
                    BoldGW[i, j] = 1.0/( 1j*w + (q*q-Para.EF) + Static[i] + (SigmaW[i,j]-dMu) )


            dG_W = BoldGW - BareGW

            dG_T, _ = Fourier.SpectralW2T(dG_W)

            fname = "green_order{0}.data".format(o)
            with open(os.path.join(selfFolder, fname), "w") as f:
                for k in range(Para.MomGridSize):
                    for t in range(Para.TauGridSize):
                        f.write("{0} ".format(dG_T[k, t].real))
                f.write("\n")
            with open(os.path.join(selfFolder, "greenList_order{0}.data".format(o)), "a+") as f:
                for k in range(Para.MomGridSize):
                    for t in range(Para.TauGridSize):
                        f.write("{0} ".format(dG_T[k, t].real))
                f.write("\n")

            if o == Para.Order:
                fname = "green.data".format(o)
                with open(os.path.join(selfFolder, fname), "w") as f:
                    for k in range(Para.MomGridSize):
                        for t in range(Para.TauGridSize):
                            f.write("{0} ".format(dG_T[k, t].real))
                    f.write("\n")
        # print("Maximum Error of \delta G: ", np.amax(abs(dG_T - dG_Tp)))

        time.sleep(2)

        flag = np.array([step/1000000 >= Para.TotalStep for step in Step])
        if np.all(flag == True):
            print("End of Simulation!")
            sys.exit(0)
