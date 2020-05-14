#!/usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import sys


def FermiKernel(w, t, beta):
    x = beta*w/2
    y = 2*t/beta-1
    if x > 100:
        return np.exp(-x*(y+1.))
    if x < -100:
        return np.exp(x*(1.0-y))
    return np.exp(-x*y)/(2*np.cosh(x))


def BoseKernel(w, t, beta):
    x = beta*w/2
    y = 2*t/beta-1
    if x > 200:
        return w*np.exp(-x*(y+1.))
    if x < -200:
        return -w*np.exp(x*(1.0-y))
    return w*np.exp(-x*y)/(2*np.sinh(x))


def BuildBasis(TauGrid, Beta, N, Type):
    Nw = 1000
    w = np.linspace(-100, 100, Nw)
    Nt = len(TauGrid)
    kMatrix = np.zeros([Nw, Nt])
    if Type == "Fermi":
        for i in range(len(w)):
            kMatrix[i, :] = FermiKernel(w[i], TauGrid, Beta)
    elif Type == "Bose":
        for i in range(len(w)):
            kMatrix[i, :] = BoseKernel(w[i], TauGrid, Beta)
    else:
        print "Not implemented!"
        sys.exit(0)

    u, s, v = linalg.svd(kMatrix)
    v_inv = linalg.inv(v)

    # Basis = []
    # for i in range(N):
    #     Basis.append(list(v_inv[:, i]))
    return v[0:N, :], v_inv[:, 0:N]


if __name__ == "__main__":
    beta = 10.0
    MaxTauBin = 1024
    Nt = MaxTauBin
    t = np.linspace(0, beta, Nt+1)
    TauGrid = np.array([e+1.0/Nt/2 for e in t[:-1]])
    svd = BuildBasis(TauGrid, beta, 15, "Fermi")
    print svd
