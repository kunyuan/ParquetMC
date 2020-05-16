import numpy as np
from utility import *
import numpy.linalg as linalg
from color import *


class grid:
    def __init__(self, para):
        self.Para = para
        self.TauSize = para.TauGridSize
        self.MomSize = para.MomGridSize
        self.AngSize = para.AngGridSize

        #### Grid ############################
        self.TauGrid = self.__TauGrid()
        self.MomGrid = self.__MomGrid()
        self.AngGrid = self.__AngleGrid()

    def __TauGrid(self):
        ########### logirthmic grid #####################
        # N = self.Para.TauGridSize/2
        # a = self.Para.EF*self.Para.Beta/6.0
        # b = a/(N-1)
        # eps = 1.0e-8/self.Para.Beta

        # TauGrid = np.zeros(self.Para.TauGridSize)

        # print a, b
        # assert a > 1, "Coeff a must be much larger than 1!"
        # assert b > 0, "Coeff b must be larger than 0!"

        # for i in range(N):
        #     TauGrid[i] = 0.5*np.exp(-a)*(np.exp(b*i)-1.0)+eps
        #     # print i, b*i, a
        #     TauGrid[2*N-1-i] = 1.0-TauGrid[i]

        # TauGrid *= self.Para.Beta

        ########## uniform grid  #########################
        arr = np.array(range(self.TauSize))
        TauGrid = arr*self.Para.Beta/(self.TauSize-1)
        TauGrid[0] += 1.0e-8
        TauGrid[-1] -= 1.0e-8
        return TauGrid

    def __MomGrid(self):
        arr = np.array(range(self.MomSize))
        KGrid = arr*self.Para.MaxExtMom/self.MomSize+1.0e-8
        return KGrid

    def __AngleGrid(self):
        arr = np.array(range(self.AngSize))
        AngGrid = arr*2.0/self.AngSize-1.0
        return AngGrid


class spectral:
    def __init__(self, para):
        self.Para = para
        ### Spectral represetation  ##########
        self.TauBasisSize = para.TauBasisSize
        self.RealFreqSize = para.RealFreqGridSize
        self.RealFreqGrid = np.linspace(
            -para.MaxRealFreq, para.MaxRealFreq, para.RealFreqGridSize)

    def TauBasis(self, TauGrid, Type, Num):
        """v[0:Num, :] maps N coefficients to Tau, v.T[:, 0:Num] maps Tau to coefficients"""
        _, s, v = self.TauKernel(TauGrid, Type)
        print yellow("Smallest Singular Value: {0} ".format(s[-1]))
        return v.T[:, :Num]

    def MatFreqBasis(self, MatFreqGrid, Type, Num):
        """v[0:Num, :] maps Tau to N coefficients, v.T[:, 0:Num] maps coefficients to Tau"""
        _, s, v = self.MatFreqKernel(MatFreqGrid, Type)
        print yellow("Smallest Singular Value: {0} ".format(s[-1]))
        return v.T[:, :Num]

    def TauKernel(self, TauGrid, Type):
        kernel = np.zeros([self.RealFreqSize, len(TauGrid)])
        for i, w in enumerate(self.RealFreqGrid):
            kernel[i, :] = self.Kernel(w, TauGrid, self.Para.Beta, Type)
        # multiply the kernel with Delta \omega
        kernel *= 2.0*self.Para.MaxRealFreq/self.RealFreqSize/2.0/np.pi
        u, s, v = linalg.svd(kernel)
        return u, s, v

    def MatFreqKernel(self, MatFreqGrid, Type):
        """v[:, 0:N] gives the first N basis for Matsubara frequency axis"""
        kernel = np.zeros([self.RealFreqSize, len(MatFreqGrid)], dtype=complex)
        for i, w in enumerate(self.RealFreqGrid):
            kernel[i, :] = 1.0/(1j*MatFreqGrid+w)
        kernel *= 2.0*self.Para.MaxRealFreq/self.RealFreqSize/2.0/np.pi
        u, s, v = linalg.svd(kernel)
        return u, s, v

    def Kernel(self, w, t, beta, Type):
        if Type == "Fermi":
            x, y = beta*w/2, 2*t/beta-1
            if x > 100:
                return np.exp(-x*(y+1.))
            elif x < -100:
                return np.exp(x*(1.0-y))
            else:
                return np.exp(-x*y)/(2*np.cosh(x))
        else:
            print "Not implemented!"
            raise ValueError


def Kernel(w, t, beta, Type):
    if Type == "Fermi":
        x, y = beta*w/2, 2*t/beta-1
        if x > 100:
            return np.exp(-x*(y+1.))
        elif x < -100:
            return np.exp(x*(1.0-y))
        else:
            return np.exp(-x*y)/(2*np.cosh(x))
    else:
        print "Not implemented!"
        raise ValueError


def TauKernel(Beta, TauGrid, RealFreqGrid, Type):
    dRealFreq = (RealFreqGrid[-1]-RealFreqGrid[0])/len(RealFreqGrid)
    kernel = np.zeros([len(RealFreqGrid), len(TauGrid)])
    for i, w in enumerate(RealFreqGrid):
        kernel[i, :] = Kernel(w, TauGrid, Beta, Type)
    # multiply the kernel with Delta \omega
    return kernel*dRealFreq/2.0/np.pi


def MatFreqKernel(Beta, MatFreqGrid, RealFreqGrid, Type):
    dRealFreq = (RealFreqGrid[-1]-RealFreqGrid[0])/len(RealFreqGrid)
    kernel = np.zeros([len(RealFreqGrid), len(MatFreqGrid)], dtype=complex)
    for i, w in enumerate(RealFreqGrid):
        kernel[i, :] = 1.0/(1j*MatFreqGrid+w)
    return kernel*dRealFreq/2.0/np.pi


def FitData(Data, Num, v):
    assert Num < v.shape[0], "Number of basis is too large!"
    coef = np.dot(Data, v.T[:, :Num])
    fitted = np.dot(coef, v[:Num, :])
    print "Max of |Sigma-Fitted Sigma|: ", np.amax(abs(Data-fitted))
    return fitted, coef


def Data2Spectral(Data, Num, u, s, v):
    ss = np.zeros((v.shape[0], u.shape[0]))
    assert Num < u.shape[0], "Number of basis is too large!"
    assert Num < v.shape[0], "Number of basis is too large!"
    for i in range(Num):
        ss[i, i] = 1.0/s[i]
    InvKernel = np.dot(np.dot(v.T, ss), u.T)
    return np.dot(Data, InvKernel)


def Spectral2Data(spectral, Num, u, s, v):
    ss = np.zeros((u.shape[0], v.shape[0]))
    assert Num < u.shape[0], "Number of basis is too large!"
    assert Num < v.shape[0], "Number of basis is too large!"
    for i in range(Num):
        ss[i, i] = s[i]
    Kernel = np.dot(np.dot(u, ss), v)
    return np.dot(spectral, Kernel)


if __name__ == "__main__":
    Para = param()
    Grid = grid(Para)
    Spectral = spectral(Para)

    with open("grid.data", "w") as f:
        f.writelines("{0} #TauGrid\n".format(Grid.TauSize))
        for t in Grid.TauGrid:
            f.write("{0} ".format(t))
        f.write("\n\n")

        f.writelines("{0} #MomGrid\n".format(Grid.MomSize))
        for k in Grid.MomGrid:
            f.write("{0} ".format(k))
        f.write("\n\n")

        f.writelines("{0} #AngGrid\n".format(Grid.AngSize))
        for a in Grid.AngGrid:
            f.write("{0} ".format(a))
        f.write("\n\n")

        basis = Spectral.TauBasis(Grid.TauGrid, "Fermi", Para.TauBasisSize)
        f.writelines("{0} #FermiTauBasis\n".format(Para.TauBasisSize))
        for t in range(Para.TauGridSize):
            for b in range(Para.TauBasisSize):
                f.write("{0} ".format(basis[t, b]))
        f.write("\n\n")

    pass
