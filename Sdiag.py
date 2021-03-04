#!/usr/bin/env python3
from utility.IO import *
from utility.plot import *
import utility.polar0 as polar
import utility.dlr_boson as dlr
import numpy as np
import grid
import scipy.linalg as slinalg
from scipy.linalg import lu_factor, lu_solve
from scipy.integrate import simps

# F+ and F-
# Fs, Fa = -0.5, -0.3

Para = param("./")

print(polar.Polarisi(0.5, 0.0, Para.EF))
print(polar.Polarisi(0.5, 0.1, Para.EF))

############# Construct Tau and Mom grids ################
T = grid.Tau()
# double beta, int _size, double scale
T.build(Para.Beta, Para.TauGridSize, 6.0/Para.EF)
# print(T.grid)
K = grid.BoseK()
# double kF, double maxK, int _size, double scale
# K.build(Para.kF, Para.MaxExtMom,
#         Para.MomGridSize, np.sqrt(1.0 / Para.Beta) * 1.5)
K.build(Para.kF/2, Para.MaxExtMom,
        Para.MomGridSize, 1.0 / Para.kF / 2.0)

# print(K.grid)

N = 10000
wnlist = np.array(range(-N, N))*2.0*np.pi/Para.Beta
w1 = np.pi/Para.Beta
weight = np.zeros(K.size)+0.0*1j

# particle-particle
for ki, k in enumerate(K.grid):
    weight[ki] = 0.0
    ek = k*k-Para.kF**2
    for n, wn in enumerate(wnlist):
        # w2 = (8.0*np.pi/(k*k+Para.Mass2-polar.Polarisi(k, wn, Para.EF)))**2
        # w2 = (8.0*np.pi/(k*k+Para.Mass2))**2
        w2 = 1.0
        weight[ki] += 1.0/(1j*(wn+w1)-ek)/(-1j*(wn+w1)-ek)*w2
    weight[ki] *= 4.0*np.pi*k*k/(2.0*np.pi)**3*Para.Beta
    print(k, weight[ki])

# bubble
# for ki, k in enumerate(K.grid):
#     weight[ki] = 0.0
#     ek = k*k-Para.kF**2
#     for n, wn in enumerate(wnlist):
#         # w2 = (8.0*np.pi/(k*k+Para.Mass2-polar.Polarisi(k, wn, Para.EF)))**2
#         weight[ki] += 1.0/(1j*wn-ek)/(1j*wn-ek)
#     weight[ki] *= 4.0*np.pi*k*k/(2.0*np.pi)**3*Para.Beta*2
#     print(k, weight[ki])

I = simps(weight, K.grid)
print("Integral: ", I*Para.Nf)

ek = np.array(K.grid)**2-Para.kF**2
w2 = (8.0*np.pi/(np.array(K.grid)**2+Para.Mass2))**2
benchmark = (1.0-np.exp(-2.0*ek*Para.Beta))/2.0 / \
    ek/(1.0+np.exp(-ek*Para.Beta))**2
benchmark *= 4.0*np.pi*k*k/(2.0*np.pi)**3
I = simps(benchmark, K.grid)
print("benchmark Integral: ", I*Para.Nf)

plt.figure()
Errorbar(K.grid[:], weight[:].real, label="freq")
Errorbar(K.grid[:], benchmark[:].real, label="benchmark")
plt.legend()
plt.show()
