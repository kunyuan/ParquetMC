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
parser.add_argument("folder1")
parser.add_argument("folder2")
args = parser.parse_args()

folderA = args.folder1
folderF = args.folder2
print("Folders : " + folderA + ", " + folderF)

legendreL = [0,1,2,3]

Para = param(folderA)
# 0: I, 1: T, 2: U, 3: S
Channel = [0, 1, 2, 3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
ChanColor = {0: "k", 1: "r", 2: "b", 3: "g"}
# 0: total, 1: order 1, ...
Order = range(Para.Order+1)
Irreducible = False

shape = (Para.Order+1, 4, Para.AngGridSize, Para.MomGridSize, 2)
Data, Norm, Step, Grid = LoadFile(folderA, "vertex_pid[0-9]+.dat", shape)

AngGrid = Grid["AngleGrid"]
MomGrid = Grid["KGrid"]
Angle = np.arccos(AngGrid)

# Data shape : (pid numbers, order, chan, AngleGrid, KGrid)

def PrintInfo(Channel, Data, DataErr):
    Data = -np.copy(Data)
    DataErr = np.copy(DataErr)

    Data *= Para.Nf
    DataErr *= Para.Nf

    # print Data.shape, DataErr.shape

    print("{0}     Q/kF,    Data,    Error".format(Channel))
    print("As: {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[0], DataErr[0]))
    print("Aa:  {0:6.2f}, {1:10.6f}, {2:10.6f}".format(
        MomGrid[0], Data[1], DataErr[1]))


def SpinMapping(Data):
    d = np.copy(Data)
    d[..., 0] += d[..., 1]/Para.Spin
    d[..., 1] /= Para.Spin
    return d


def Bare(angle, Lambda):

    Bare = np.zeros([Para.AngGridSize, 2])
    if Irreducible == False:
        Bare[:, 0] += -8.0*np.pi/(Para.Mass2+Lambda)*Para.Nf
    Bare[:, 1] += +8.0 * np.pi / \
        ((2.0*Para.kF*np.sin(angle/2.0))**2+Para.Mass2+Lambda)*Para.Nf
    Bare = SpinMapping(Bare)
    return Bare


bare = Bare(Angle, 0.0)

Data = [np.sum(d[1:Para.Order+1, ...], axis=0) for d in Data]

#---------------- 
Adata = [SpinMapping(np.sum(d[:, :, 0, :], axis=0))*Para.Nf for d in Data]
avg, err = Estimate(Adata, Norm)
bareLambda = Bare(Angle, Para.Lambda)

Adata_s = -(avg[:, 0]+bareLambda[:, 0])
Adata_a = -(avg[:, 1]+bareLambda[:, 1])
Aerr_s = -err[:, 0]
Aerr_a = -err[:, 1]

Bl = LegendreCoeff(Adata_s, AngGrid, legendreL)
Cl = LegendreCoeff(Adata_a, AngGrid, legendreL)
BlErr = LegendreCoeff(Aerr_s, AngGrid, legendreL)
ClErr = LegendreCoeff(Aerr_a, AngGrid, legendreL)


#----------------------------------------------
Para = param(folderF)
Data, Norm, Step, Grid = LoadFile(folderF, "vertex_pid[0-9]+.dat", shape)

AngGrid = Grid["AngleGrid"]
MomGrid = Grid["KGrid"]
Angle = np.arccos(AngGrid)

Data = [np.sum(d[1:Para.Order+1, ...], axis=0) for d in Data]

#---------------- 
Fdata = [SpinMapping(np.sum(d[:, :, 0, :], axis=0))*Para.Nf for d in Data]
avg, err = Estimate(Fdata, Norm)
Irreducible = True
bareLambda = Bare(Angle, Para.Lambda)
Fdata_s = -(avg[:, 0]+bareLambda[:, 0])
Fdata_a = -(avg[:, 1]+bareLambda[:, 1])
Ferr_s = -err[:, 0]
Ferr_a = -err[:, 1]

Fl = LegendreCoeff(Fdata_s, AngGrid, legendreL)
Zl = LegendreCoeff(Fdata_a, AngGrid, legendreL)
FlErr = LegendreCoeff(Ferr_s, AngGrid, legendreL)
ZlErr = LegendreCoeff(Ferr_a, AngGrid, legendreL)

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
        print("order-" + str(k) + ":  ", Bl[k], "+-", abs(Errl[k]) )


print("Fl: ")
AFprint(Fl, FlErr)
print("\nZl:")
AFprint(Zl, ZlErr)

print("\n\nrenormalized Fl :")
AFprint(renormalize(Fl), renormalize_error(Fl, FlErr))
print("Bl:")
AFprint(Bl, BlErr)

print("\nrenormalized Zl:")
AFprint(renormalize(Zl), renormalize_error(Zl, ZlErr))
print("Cl:")
AFprint(Cl, ClErr)
