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
        # Size = GridSize-1
        # arr = np.array(range(Size+1))
        # TauGrid = arr*1.0
        # c = Para.Beta*Para.EF/3.0
        # for i in range(Size/2):
        #     TauGrid[i] = Para.Beta/2.0 * \
        #         np.exp(c*(float(i)/Size-0.5)) - \
        #         Para.Beta/2.0*np.exp(-c*0.5)+1.0e-8
        # for i in range(Size/2, Size+1):
        #     TauGrid[i] = Para.Beta-Para.Beta/2.0 * \
        #         np.exp(c*(0.5-float(i)/Size))+Para.Beta/2.0 * \
        #         np.exp(-c*0.5)-1.0e-8

        ########## uniform grid  #########################
        arr = np.array(range(self.TauSize))
        TauGrid = arr*self.Para.Beta/(self.TauSize-1)
        TauGrid[0] += 1.0e-8
        TauGrid[-1] -= 1.0e-8
        return TauGrid

    def __MomGrid(self):
        arr = np.array(range(self.MomSize))
        KGrid = arr*Para.MaxExtMom/self.MomSize+1.0e-8
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
        """u.T[0:Num, :] gives the first N basis for Tau axis"""
        u, s, _ = self.TauKernel(TauGrid, Type)
        print yellow("Singular Value (0-{0}): ".format(Num))
        print s[:Num]
        return u.T[:Num, :]

    def MatFreqBasis(self, MatFreqGrid, Type, Num):
        """u.T[0:Num, :] gives the first N basis for Tau axis"""
        u, s, _ = self.MatFreqKernel(MatFreqGrid, Type)
        print yellow("Singular Value (0-{0}): ".format(Num))
        print s[:Num]
        return u.T[:Num, :]

    def TauKernel(self, TauGrid, Type):
        kernel = np.zeros([len(TauGrid), self.RealFreqSize])
        for i, w in enumerate(self.RealFreqGrid):
            kernel[:, i] = self.Kernel(w, TauGrid, Para.Beta, Type)
        # multiply the kernel with Delta \omega
        kernel *= 2.0*self.Para.MaxRealFreq/self.RealFreqSize
        u, s, v = linalg.svd(kernel)
        return u, s, v

    def MatFreqKernel(self, MatFreqGrid, Type):
        """v[:, 0:N] gives the first N basis for Matsubara frequency axis"""
        kernel = np.zeros([len(MatFreqGrid), self.RealFreqSize], dtype=complex)
        for i, w in enumerate(self.RealFreqGrid):
            kernel[:, i] = 1.0/(1j*MatFreqGrid-w)
        kernel *= 2.0*self.Para.MaxRealFreq/self.RealFreqSize
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
        for b in range(Para.TauBasisSize):
            for t in range(Para.TauGridSize):
                f.write("{0} ".format(basis[b, t]))
        f.write("\n\n")

    pass
