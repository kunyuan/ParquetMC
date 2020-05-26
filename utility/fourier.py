import numpy as np
from scipy import integrate
from scipy import interpolate
import spectral
import numpy.linalg as linalg


class fourier:
    def __init__(self, TauGrid, WnGrid, Beta):
        self.TauGrid = TauGrid
        self.WnGrid = WnGrid
        self.TauSize = len(TauGrid)
        self.WnSize = len(WnGrid)
        self.Beta = Beta

        self.uTauGrid = np.linspace(0.0, Beta, 4096)  # uniform grid
        print "TauGrid", self.TauGrid[0], Beta-self.TauGrid[-1]
        self.uTauGrid[0] = 2.e-8
        self.uTauGrid[-1] = Beta-4.e-7

    def InitializeKernel(self, MaxRealFreq, RealFreqSize, Type, Threshold):
        self.RealFreqGrid = np.linspace(-MaxRealFreq,
                                        MaxRealFreq, RealFreqSize)
        self.RealFreqSize = RealFreqSize
        self.Type = Type
        self.KernelT = spectral.TauKernel(
            self.Beta, self.TauGrid, self.RealFreqGrid, self.Type)

        self.uKernelT = spectral.TauKernel(
            self.Beta, self.uTauGrid, self.RealFreqGrid, self.Type)
        self.uInvKernelT = linalg.pinv(self.uKernelT, Threshold)

        self.KernelW = spectral.MatFreqKernel(
            self.Beta, self.WnGrid, self.RealFreqGrid, self.Type)
        self.InvKernelW = linalg.pinv(self.KernelW, Threshold)

    def SpectralT2W(self, dataT):

        f = interpolate.interp1d(self.TauGrid, dataT, axis=-1)
        dataT = f(self.uTauGrid)

        # spectral = np.dot(dataT, self.__pinv(ut, st, vt, w, threshold))
        spectral = np.dot(dataT, self.uInvKernelT)
        dataW = np.dot(spectral, self.KernelW)
        # kernel = linalg.pinv(linalg.pinv(
        #     self.KernelW, 1.0e-14), 1.0e-14)
        # dataW = np.dot(spectral, kernel)
        return dataW, spectral

    def SpectralW2T(self, dataW):
        spectral = np.dot(dataW, self.InvKernelW)
        print self.InvKernelW.shape
        print dataW.shape
        dataT = np.dot(spectral, self.uKernelT)
        f = interpolate.interp1d(self.uTauGrid, dataT, axis=-1)
        sTauGrid = np.array(self.TauGrid)
        sTauGrid[0] = 3.0e-8
        sTauGrid[-1] = self.Beta-5.0e-7
        # print self.uTauGrid[0], self.uTauGrid[-1]
        # print self.TauGrid[0], self.TauGrid[-1]
        dataT = f(sTauGrid)
        return dataT, spectral

    # def __pinv(self, u, s, v, w, threshold):
    #     ss = np.zeros([v.shape[0], u.shape[0]])
    #     print np.dot(v, v.T)
    #     # print np.dot(np.dot(v, w), v.T)
    #     for i in range(np.min([ss.shape[0], ss.shape[1]])):
    #         if s[i] > threshold:
    #             ss[i, i] = 1.0/s[i]/w[i, i]
    #     # return np.dot(np.dot(v.T, ss), u.T)
    #     return np.dot(np.dot(np.dot(w, v.T), ss), u.T)

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

        dataT = np.zeros(Tshape, dtype=complex)
        for i, tau in enumerate(self.TauGrid):
            phase = np.exp(1j*self.WnGrid*tau)
            dt = dataW[..., :]*phase
            dataT[..., i] = np.sum(dt, axis=-1)
        return dataT/self.Beta
