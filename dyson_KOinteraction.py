#!/usr/bin/env python3
from utility.IO import *
from utility.plot import *
import utility.polar0 as polar
import utility.dlr_boson as dlr
import numpy as np
import grid
import scipy.linalg as slinalg
from scipy.linalg import lu_factor, lu_solve

# F+ and F-
Fs, Fa = -0.5, -0.3

Para = param("./")

polar.Polarisi(1.e-6, 0.0, Para.EF)

############# Construct Tau and Mom grids ################
T = grid.Tau()
# double beta, int _size, double scale
T.build(Para.Beta, Para.TauGridSize, 6.0/Para.EF)
# print(T.grid)
Kbose = grid.BoseK()
# double kF, double maxK, int _size, double scale
Kbose.build(Para.kF, Para.MaxExtMom,
            Para.MomGridSize, np.sqrt(1.0/Para.kF))
print(Kbose.grid)

############  DLR basis   #################################
Wmax = 6.0
Nw = 1024
dW = Wmax/Nw
Nt = int(Para.Beta*Wmax*100)
Nwn = int(Wmax*Para.Beta/2.0/np.pi*100)
eps = 1.0e-12
# print("Nt: ", Nt)
# print("Nwn: ", Nwn)
k, wGrid, wnGrid, tGrid = dlr.getDLR(Para.Beta, Wmax, Nw, Nwn, Nt, eps)
for i in range(k):
    print(i, wGrid[i], tGrid[i], wnGrid[i])

PolarW = np.zeros([Kbose.size, k])
wnlist = wnGrid*2.0*np.pi/Para.Beta
for qi, q in enumerate(Kbose.grid):
    PolarW[qi, :] = polar.Polarisi(q, wnlist, Para.EF)

coeff = np.zeros([Kbose.size, k])
PolarBench = np.zeros([Kbose.size, k])
KerW = dlr.getKerW(wnGrid, wGrid, Para.Beta)
lu, piv = lu_factor(KerW)

# calculate the dlr coefficients
for qi in range(Kbose.size):
    coeff[qi, :] = lu_solve((lu, piv), PolarW[qi, :])
    PolarBench[qi, :] = KerW @ coeff[qi, :]

print("Maximum fitting error: ", np.max(abs(PolarBench-PolarW)))

############ Plot Polarization in Matsubara frequency ################
# plt.figure()
# for wi, w in enumerate(wnlist):
#     Errorbar(Kbose.grid, Polar[:, wi], label=f"{wnGrid[wi]}")
#     # Errorbar(Kbose.grid, PolarBench[:, wi], label=f"{wnGrid[wi]}")
# plt.legend()
# plt.show()

# calculate the polarization in tau
KerT = dlr.getKerT(T.grid, wGrid, Para.Beta)
PolarT = np.zeros([Kbose.size, T.size])
for qi, q in enumerate(Kbose.grid):
    PolarT[qi, :] = KerT @ coeff[qi, :]

print("Maximu error in Tau:", np.max(abs(PolarW[0, 0]-PolarT[0, :]*Para.Beta)))

############ Plot Polarization in Tau ################
# plt.figure()
# for qi, q in enumerate(Kbose.grid[::2]):
#     Errorbar(T.grid, PolarT[qi, :], label=f"{q}")
# plt.legend()
# plt.show()

########### KO interaction ###########################


# def Gp(q): return gp*(q/Para.qTF)**2
# def Gm(q): return gm*(q/Para.qTF)**2


# KO interaction:
# Rs=(v-fs)/(1-(v-fs)*Pi0)
# Ra=-fa/(1+fa*Pi0)
def fs(q, Fs): return Fs/Para.Nf
def fa(q, Fa): return Fa/Para.Nf


# plot the denorminator of Rs
# Fs = np.array([0.6, 0.8, 1.0, 1.5, 1.75, 2.0])

# denorm = np.zeros([len(Fs), Kbose.size])
# for qi, q in enumerate(Kbose.grid):
#     inv = q*q/(8.0*np.pi)
#     denorm[:, qi] = (inv-(1.0-inv*fs(q, Fs))*PolarW[qi, 0])/abs(PolarW[0, 0])

# plt.figure()
# for fi, f in enumerate(Fs):
#     Errorbar(Kbose.grid, denorm[fi, :], label=f"{f}")
# plt.grid()
# plt.ylim([0.0, 2.0])
# plt.xlim([0.0, 2.0])
# plt.legend()
# plt.show()


# we only calculate the part which is continuous in tau
dRsw = np.zeros([Kbose.size, len(wnGrid)])
dRaw = np.zeros([Kbose.size, len(wnGrid)])
dRsT = np.zeros([Kbose.size, T.size])
dRaT = np.zeros([Kbose.size, T.size])


for qi, q in enumerate(Kbose.grid):
    inv = (q*q+0.01)/(8.0*np.pi)  # add a small mass to regularize everything
    # dRs=(v+fs)^2*Pi0/(1-(v+fs)*Pi0)
    denorm = inv-(1.0+inv*fs(q, Fs))*PolarW[qi, :]
    dRsw[qi, :] = (1.0+inv*fs(q, Fs))**2*PolarW[qi, :]/denorm/inv
    coeff[qi, :] = lu_solve((lu, piv), dRsw[qi, :])
    dRsT[qi, :] = KerT @ coeff[qi, :]

    # dRa=fa^2*Pi0/(1-fa*Pi0)
    denorm = 1.0-fa(q, Fa)*PolarW[qi, :]
    dRaw[qi, :] = fa(q, Fa)**2*PolarW[qi, :]/denorm
    coeff[qi, :] = lu_solve((lu, piv), dRaw[qi, :])
    dRaT[qi, :] = KerT @ coeff[qi, :]

# Ra

########### Plot Polarization in Tau ################
plt.figure()
for qi, q in enumerate(Kbose.grid[::6]):
    Errorbar(T.grid, dRsT[qi, :], label=f"{q}")
plt.title("dRs")
plt.legend()
plt.show()

plt.figure()
for qi, q in enumerate(Kbose.grid[::6]):
    Errorbar(T.grid, dRaT[qi, :], label=f"{q}")
plt.title("dRa")
plt.legend()
plt.show()
