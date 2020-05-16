import numpy as np
import grid
from scipy import integrate
from scipy.interpolate import UnivariateSpline
import numpy.linalg as linalg


class fourier:
    def __init__(self, TauGrid, WnGrid, Beta):
        self.TauGrid = TauGrid
        self.WnGrid = WnGrid
        self.TauSize = len(TauGrid)
        self.WnSize = len(WnGrid)
        self.Beta = Beta

    def InitializeKernel(self, MaxRealFreq, RealFreqSize, Type):
        self.RealFreqGrid = np.linspace(-MaxRealFreq,
                                        MaxRealFreq, RealFreqSize)
        self.RealFreqSize = RealFreqSize
        self.Type = Type
        self.KernelT = grid.TauKernel(
            self.Beta, self.TauGrid, self.RealFreqGrid, self.Type)
        self.usvT = linalg.svd(self.KernelT)

        self.KernelW = grid.MatFreqKernel(
            self.Beta, self.WnGrid, self.RealFreqGrid, self.Type)
        self.usvW = linalg.svd(self.KernelW)

    def SpectralT2W(self, dataT, weight=None, threshold=0.0):
        ut, st, vt = self.usvT

        if weight == None:
            w = np.identity(self.TauSize)
        else:
            assert len(weight) == self.TauSize, "There must be {0} number of weights!".format(
                self.TauSize)
            w = np.diag(weight)

        spectral = np.dot(dataT, self.__pinv(ut, st, vt, w, threshold))
        # spectral = np.dot(dataT, linalg.pinv(self.KernelT, threshold))
        dataW = np.dot(spectral, self.KernelW)
        return dataW, spectral

    def __pinv(self, u, s, v, w, threshold):
        ss = np.zeros([v.shape[0], u.shape[0]])
        for i in range(np.min([ss.shape[0], ss.shape[1]])):
            if s[i] > threshold:
                ss[i, i] = 1.0/s[i]
        return np.dot(np.dot(v.T, ss), u.T)
        # return np.dot(np.dot(np.dot(w, v.T), ss), u.T)

    def naiveT2W(self, dataT):
        """do fourier transformation on the last axis"""
        assert dataT.shape[-1] == self.TauSize, "The Tau axis must have {0} elements.".format(
            self.TauSize)

        Wshape = np.array(dataT.shape)
        Wshape[-1] = self.WnSize

        dataW = np.zeros(Wshape, dtype=complex)
        for i, freq in enumerate(self.WnGrid):
            phase = np.exp(-1j*self.TauGrid*freq)
            dw = dataT[..., :]*phase
            dataW[..., i] = integrate.trapz(dw[..., :], self.TauGrid)
            # SigmaW[i] = integrate.simps(dw, self.TauGrid)
        return dataW

    def naiveW2T(self, dataW):
        """do fourier transformation on the last axis"""
        assert dataW.shape[-1] == self.WnSize, "The Tau axis must have {0} elements.".format(
            self.WnSize)

        Tshape = np.array(dataW.shape)
        Tshape[-1] = self.TauSize

        dataT = np.zeros(Tshape)
        for i, tau in enumerate(self.TauGrid):
            phase = np.exp(1j*self.WnGrid*tau)
            dt = dataW[..., :]*phase
            dataT[..., i] = np.sum(dt, axis=-1)
        return dataT/self.Beta
