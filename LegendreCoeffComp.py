#!/usr/bin/env python3
from utility.IO import *
from utility.angle import LegendreCoeff
import numpy as np
import sys
import argparse
# import matplotlib.pyplot as plt
# import matplotlib as mat
# mat.rcParams.update({'font.size': 16})
# mat.rcParams["font.family"] = "Times New Roman"
# size = 12

parser = argparse.ArgumentParser("Specify some parameters.")
parser.add_argument("--folderA", "-A", required=True)
parser.add_argument("--folderF", "-F", required=True)
args = parser.parse_args()

folderA = args.folderA
folderF = args.folderF
print("Folders : " + folderA + ", ")


Type = "As_Aa"
# Type = "Dir_Ex"
# Type = "Gamma_4spin"
legendreL = [0]

Para = param(folderA)

# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
ChanColor = {0: "k", 1: "r", 2: "b", 3: "g"}
# 0: total, 1: order 1, ...
Order = range(Para.Order+1)
Irreducible = True

shape = (Para.Order+1, 4, Para.AngGridSize, Para.MomGridSize, 2)
Data, Norm, Step, Grid = LoadFile(folderA, "vertex_pid[0-9]+.dat", shape)

AngGrid = Grid["AngleGrid"]         # cos
MomGrid = Grid["KGrid"]
Angle = np.arccos(AngGrid)

# Data shape : (pid numbers, order, chan, AngleGrid, KGrid, 2)


def SpinMapping(Data):
    d = np.copy(Data)
    if Type == "As_Aa":
        d[..., 0] += d[..., 1]/Para.Spin
        d[..., 1] /= Para.Spin
    elif Type == "Gamma_4spin":
        e = d[..., 0] + d[..., 1]
        d[..., 1] = d[..., 0]
        d[..., 0] = e
    elif Type == "Dir_Ex":
        pass
    return d


def Bare(angle, Lambda):
    Bare = np.zeros([Para.AngGridSize, 2])
    if Irreducible == False:
        Bare[:, 0] += -8.0*np.pi/(Para.Mass2+Lambda)
    Bare[:, 1] += +8.0 * np.pi / \
        ((2.0*Para.kF*np.sin(angle/2.0))**2+Para.Mass2+Lambda)
    Bare = SpinMapping(Bare)
    return Bare



def GetCoeff():
    m = 1.0/2.0
    Zname = "ZandM_order{0}.data".format(Para.Order-1)
    try:
        Z, mStar = np.loadtxt(os.path.join(folderF, "selfconsistent", Zname))
    except Exception as e:
        print("Can not load Z and m*")
        Z, mStar = 1, m
    # Z, mStar = 0.873, 0.955*0.5
    coeff = Z*Z * mStar * Para.kF / ( np.pi*np.pi)
    return coeff


bare = Bare(Angle, 0.0)

# Data shape : (pid numbers, order, chan, AngleGrid, KGrid, 2)
Data1 = [np.sum(d[1:Para.Order+1, ...], axis=0) for d in Data]
# Data shape : (pid numbers, chan, AngleGrid, KGrid, 2)


#=========================== A ===============================
Adata = [SpinMapping(np.sum(d[:, :, 0, :], axis=0)) for d in Data1]
# Adata shape : (pid numbers, AngleGrid, 2), channels are summed, momentum is set as p=0. 
avg, err = Estimate(Adata, Norm)
bareLambda = Bare(Angle, Para.Lambda)

Adata_s = -(avg[:, 0]+bareLambda[:, 0])
Adata_a = -(avg[:, 1]+bareLambda[:, 1])


# A Error 
Data2 = [np.sum(d[1:Para.Order, ...], axis=0) for d in Data]
Adata = [SpinMapping(np.sum(d[:, :, 0, :], axis=0)) for d in Data2]
avg, err = Estimate(Adata, Norm)
Adata_sErr = -(avg[:, 0]+bareLambda[:, 0])
Adata_aErr = -(avg[:, 1]+bareLambda[:, 1])


#=========================== F ===============================
Para = param(folderF)
shape = (Para.Order+1, 4, Para.AngGridSize, Para.MomGridSize, 2)
Data, Norm, Step, Grid = LoadFile(folderF, "vertex_pid[0-9]+.dat", shape)

AngGrid = Grid["AngleGrid"]
MomGrid = Grid["KGrid"]
Angle = np.arccos(AngGrid)

Data1 = [np.sum(d[1:Para.Order+1, ...], axis=0) for d in Data]
Fdata = [SpinMapping(np.sum(d[:, :, 0, :], axis=0)) for d in Data1]
avg, err = Estimate(Fdata, Norm)
bareLambda = Bare(Angle, Para.Lambda)
Fdata_s = -(avg[:, 0]+bareLambda[:, 0])
Fdata_a = -(avg[:, 1]+bareLambda[:, 1])

# F Error
Data2 = [np.sum(d[1:Para.Order, ...], axis=0) for d in Data]
Fdata = [SpinMapping(np.sum(d[:, :, 0, :], axis=0)) for d in Data2]
avg, err = Estimate(Fdata, Norm)
bareLambda = Bare(Angle, Para.Lambda)
Fdata_sErr = -(avg[:, 0]+bareLambda[:, 0])
Fdata_aErr = -(avg[:, 1]+bareLambda[:, 1])


# ============================================================================

coeff = GetCoeff()

As = LegendreCoeff(Adata_s*coeff, AngGrid, legendreL)
Aa = LegendreCoeff(Adata_a*coeff, AngGrid, legendreL)
AsErr = LegendreCoeff((Adata_s-Adata_sErr)*coeff, AngGrid, legendreL)
AaErr = LegendreCoeff((Adata_a-Adata_aErr)*coeff, AngGrid, legendreL)


Fs = LegendreCoeff(Fdata_s*coeff, AngGrid, legendreL)
Fa = LegendreCoeff(Fdata_a*coeff, AngGrid, legendreL)
# FsErr = LegendreCoeff(Ferr_s*coeff, AngGrid, legendreL)
# FaErr = LegendreCoeff(Ferr_a*coeff, AngGrid, legendreL)
FsErr = LegendreCoeff((Fdata_s-Fdata_sErr)*coeff, AngGrid, legendreL)
FaErr = LegendreCoeff((Fdata_a-Fdata_aErr)*coeff, AngGrid, legendreL)



# ============================================================================

def renormalize(ldict):
    newl = {}
    for l in ldict:
        newl[l] = ldict[l]/(1 + ldict[l]/(2*l+1))
    return newl
def renormalize_error(ldict, lerrdict):
    newl = {}
    for l in ldict:
        x = ldict[l]/(2*l+1)
        partial_Fl = 1.0/(1+x) - x/((1+x)**2)
        newl[l] = partial_Fl * lerrdict[l]
    return newl
def AFprint(Bl, Errl):
    for k in Bl.keys():
        print("legendreL l=" + str(k) + ":  ", Bl[k], "+-", abs(Errl[k]) )

print("Fs: ")
AFprint(Fs, FsErr)
print("Fa:")
AFprint(Fa, FaErr)

print("----------------------------------\nAs:")
AFprint(As, AsErr)
print("Aa:")
AFprint(Aa, AaErr)

print("----------------------------------\nRenoemalized from Fs:")
AFprint(renormalize(Fs), renormalize_error(Fs, FsErr))
print("Renoemalized from Fa:")
AFprint(renormalize(Fa), renormalize_error(Fa, FaErr))



# print("\n\nrenormalized Fs :")
# AFprint(renormalize(Fs), renormalize_error(Fs, FsErr))
# print("As:")
# AFprint(As, AsErr)

# print("\nrenormalized Fa:")
# AFprint(renormalize(Fa), renormalize_error(Fa, FaErr))
# print("Aa:")
# AFprint(Aa, AaErr)
