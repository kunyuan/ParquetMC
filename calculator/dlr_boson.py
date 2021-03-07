import numpy as np
import scipy.linalg as linalg
import scipy.linalg.interpolative as inter


def kernelT(t, w, beta):
    # print(t)
    return np.exp(-w*t)+np.exp(-w*(beta-t))


def getKerT(tlist, wlist, beta):
    tlist = np.array(tlist)
    wlist = np.array(wlist)
    KerT = np.zeros([len(tlist), len(wlist)])
    # for ti, t in enumerate(tlist):
    for wi, w in enumerate(wlist):
        KerT[:, wi] = kernelT(tlist, w, beta)
    return KerT


def kernelW(nlist, w, beta):
    wn = 2.0*np.pi*nlist/beta
    x = w*beta
    if x < 1.0e-5:
        kernel = np.zeros(len(nlist))
        for ni, flag in enumerate((nlist == 0)):
            if flag:
                kernel[ni] = (2.0-x)*beta
            else:
                kernel[ni] = 2.0*w*(x-x*x/2.0)/(wn[ni]**2+w**2)
            return kernel
    else:
        return 2.0*w*(1.0-np.exp(-x))/(wn**2+w**2)


def getKerW(nlist, wlist, beta):
    nlist = np.array(nlist)
    wlist = np.array(wlist)
    KerW = np.zeros([len(nlist), len(wlist)])
    # for ni, n in enumerate(nlist):
    for wi, w in enumerate(wlist):
        KerW[:, wi] = kernelW(nlist, w, beta)

    return KerW


def getDLR(beta, Wmax, Nw, Nwn, Nt, eps):
    dW = Wmax/Nw
    nlist = range(Nwn)
    wlist = [wi*dW for wi in range(Nw)]
    KerW = getKerW(nlist, wlist, beta)

    k, idx, _ = inter.interp_decomp(KerW, eps)
    wGrid = np.sort(idx[:k])
    wGrid = wGrid*dW

    iWn, _ = inter.interp_decomp(KerW[:, idx[:k]].T, k)
    wnGrid = np.sort(iWn[:k])

    dt = beta/(Nt-1)
    tlist = [ti*dt for ti in range(Nt)]
    KerT = getKerT(tlist, wGrid*dW, beta)
    it, _ = inter.interp_decomp(KerT.T, k)
    tGrid = np.sort(it[:k])
    tGrid = tGrid*dt
    tGrid = (tGrid+beta-tGrid[::-1])/2.0
    return k, wGrid, wnGrid, tGrid


if __name__ == "__main__":
    Wmax = 6.0
    Nw = 1024
    dW = Wmax/Nw
    beta = 25.0
    Nt = int(beta*Wmax*100)
    Nwn = int(Wmax*beta/2.0/np.pi*100)
    eps = 1.0e-16
    print("Nt: ", Nt)
    print("Nwn: ", Nwn)
    k, wGrid, wnGrid, tGrid = getDLR(beta, Wmax, Nw, Nwn, Nt, eps)
    for i in range(k):
        print(i, wGrid[i], tGrid[i], wnGrid[i])
