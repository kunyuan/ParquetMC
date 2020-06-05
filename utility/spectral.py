import numpy as np
from utility.IO import *
import numpy.linalg as linalg


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
        print(yellow(f"Smallest Singular Value: {s[-1]} "))
        return v.T[:, :Num]

    def MatFreqBasis(self, MatFreqGrid, Type, Num):
        """v[0:Num, :] maps Tau to N coefficients, v.T[:, 0:Num] maps coefficients to Tau"""
        _, s, v = self.MatFreqKernel(MatFreqGrid, Type)
        print(yellow(f"Smallest Singular Value: {s[-1]} "))
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
            print("Not implemented!")
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
    elif Type == "Bose":
        # print(w)
        # assert t > 0.0, "Tau must be positive!"
        # assert w > 0.0, "Omega must be positive!"
        return np.exp(-w*t)+np.exp(-w*(beta-t))
    else:
        print("Not implemented!")
        raise ValueError


def TauKernel(Beta, TauGrid, RealFreqGrid, Type):
    dRealFreq = (RealFreqGrid[-1]-RealFreqGrid[0])/len(RealFreqGrid)
    kernel = np.zeros([len(RealFreqGrid), len(TauGrid)])
    for i, w in enumerate(RealFreqGrid):
        kernel[i, :] = Kernel(w, TauGrid, Beta, Type)
        # print(i, w,  kernel[i, 0])
    # multiply the kernel with Delta \omega
    return kernel*dRealFreq/2.0/np.pi


def MatFreqKernel(Beta, MatFreqGrid, RealFreqGrid, Type):
    if Type == "Fermi":
        dRealFreq = (RealFreqGrid[-1]-RealFreqGrid[0])/len(RealFreqGrid)
        kernel = np.zeros([len(RealFreqGrid), len(MatFreqGrid)], dtype=complex)
        for i, w in enumerate(RealFreqGrid):
            kernel[i, :] = 1.0/(1j*MatFreqGrid+w)
        return kernel*dRealFreq/2.0/np.pi
    elif Type == "Bose":
        # require w in [0.0, +inf)
        dRealFreq = (RealFreqGrid[-1]-RealFreqGrid[0])/len(RealFreqGrid)
        kernel = np.zeros([len(RealFreqGrid), len(MatFreqGrid)], dtype=complex)
        for i, w in enumerate(RealFreqGrid):
            for j, wn in enumerate(MatFreqGrid):
                if abs(w) < 1.0e-8 and wn == 0:
                    kernel[i, j] = 2.0*Beta
                else:
                    kernel[i, j] = 2.0*w/(wn**2+w**2)*(1.0-np.exp(-Beta*w))
        return kernel*dRealFreq/2.0/np.pi
    else:
        print("Not implemented!")
        raise ValueError


def FitData(Data, Num, v):
    assert Num < v.shape[0], "Number of basis is too large!"
    coef = np.dot(Data, v.T[:, :Num])
    fitted = np.dot(coef, v[:Num, :])
    print("Max of |Sigma-Fitted Sigma|: ", np.amax(abs(Data-fitted)))
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
    Spectral = spectral(Para)

    pass
