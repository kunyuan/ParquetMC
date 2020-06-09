import numpy as np
import math


class Coeff:
    def __init__(self, b, i, l, dense2sparse):
        """idx[0] is at bound[0], idx[1] is at bound[1]"""
        self.bound = b
        self.idx = i
        self.λ = l
        self.dense2sparse = dense2sparse

        if dense2sparse == False:
            self.bound = self.bound[1], self.bound[0]
            self.idx = self.idx[1], self.idx[0]
            self.λ = -self.λ

        _l0, _l1 = 1.0, np.exp(self.λ*(self.idx[1]-self.idx[0]))
        self.b = (self.bound[1]-self.bound[0])/(_l1-_l0)
        self.a = (self.bound[0]*_l1-self.bound[1]*_l0)/(_l1-_l0)

    def floor(self, x):
        pos = self.idx[0]+1.0/self.λ*np.log((x-self.a)/self.b)
        return int(math.floor(pos))

    def grid(self, idx):
        return self.a+self.b*np.exp(self.λ*(idx-self.idx[0]))


def logGrid(segments):
    grid = []
    for coeff, rage in segments:
        for idx in rage:
            grid.append(coeff.grid(idx))
    return np.array(grid)


def tauGrid(β, Ef, size):
    assert size % 2 == 0, "size of the TauGrid must be an even number"
    λ = β/Ef/size/3.0

    c1 = Coeff([0.0, β/2.0], [0.0, size/2-0.5], λ, True)
    range1 = range(0, int(size/2))

    c2 = Coeff([β/2, β], [size/2-0.5, size-1], λ, False)
    range2 = range(int(size/2), size)

    tau = logGrid([(c1, range1), (c2, range2)])
    tau[0] = 1.0e-8
    tau[-1] = β-1.0e-8
    assert len(tau) == size, f"TauGrid size {len(tau)} is not {size}!"
    return tau


def kFermiGrid(β, Kf, MaxK, size):
    assert MaxK > Kf, "MaxK must be larger than Kf"
    kFi = size/2
    Ef = Kf**2
    λ = np.sqrt(Ef*β)/Kf

    c1 = Coeff([0.0, Kf], [0.0, kFi], λ, False)
    r1 = range(0, kFi)

    c2 = Coeff([])
