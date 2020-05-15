import os
import sys
import re
import glob
import numpy as np
from color import *

import matplotlib.pyplot as plt
import matplotlib as mat
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

##################### Global variables  #############################
Dim = 3
SpinIndex = 2
#####################################################################


def GetLine(file):
    while True:
        line = file.readline().strip()
        if len(line) > 0 and line[0] != "#":
            return line


class param:
    # Order, Beta, Rs, Mass2, Lambda, Charge2, TotalStep = [None, ]*7
    # kF, Nf, EF, Bubble = [0.0, ]*4
    InputFile = "parameter"

    def __init__(self):

        with open(self.InputFile, "r") as file:
            file.readline()  # comment line
            para = GetLine(file).split(",")
            self.Order = int(para[0])
            self.Beta = float(para[1])
            self.Rs = float(para[2])
            self.Mass2 = float(para[3])
            self.Lambda = float(para[4])
            self.Charge2 = float(para[5])
            self.TotalStep = int(para[6])

            grid = GetLine(file).split(",")
            self.TauGridSize = int(grid[0])
            self.MomGridSize = int(grid[1])
            self.AngGridSize = int(grid[2])
            self.RealFreqGridSize = int(grid[3])
            self.MaxRealFreq = float(grid[4])
            self.MaxExtMom = float(grid[5])
            self.TauBasisSize = int(grid[6])

            timer = GetLine(file).split(",")
            self.PrintTimer = int(timer[0])
            self.SaveTimer = int(timer[1])
            self.ReWeightTimer = int(timer[2])
            self.MessageTimer = int(timer[3])
            self.CollectionTimer = int(timer[4])

        if Dim == 3:
            self.kF = (9.0*np.pi/4.0)**(1.0/3.0)/self.Rs
            self.Nf = self.kF/4.0/np.pi**2*SpinIndex
        elif Dim == 2:
            self.kF = np.sqrt(2.0)/self.Rs  # 2D
            print "Not Implemented for Dimension {0}".format(Dim)
            sys.exit(0)
        else:
            print "Not Implemented for Dimension {0}".format(Dim)
            sys.exit(0)

        self.EF = self.kF**2
        self.Beta /= self.EF
        self.MaxExtMom *= self.kF
        self.MaxRealFreq *= self.EF

        print yellow("Parameters:")
        print "Rs={0}, kF={1}, EF={2}, Beta={3}, Mass2={4}, Lambda={5}\n".format(
            self.Rs, self.kF, self.EF, self.Beta, self.Mass2, self.Lambda)

        print yellow("Grid Information:")
        print "TauSize={0}, MomSize={1}, AngleSize={2}, RealFreqSize={3}, TauBasisNum={4}".format(
            self.TauGridSize, self.MomGridSize, self.AngGridSize, self.RealFreqGridSize, self.TauBasisSize)
        print "MaxRealFreq={0}, MaxExtMom={1}\n".format(
            self.MaxRealFreq, self.MaxExtMom)

        print yellow("Timer Information:")
        print "Print={0}, Save={1}, ReWeight={2}, Message={3}, Collection={4}".format(
            self.PrintTimer, self.SaveTimer, self.ReWeightTimer, self.MessageTimer, self.CollectionTimer)

# For the given path, get the List of all files in the directory tree


def getListOfFiles(dirName):
    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(dirName):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    return listOfFiles


def Estimate(Data, Weights):
    """ Return Mean and Error  with given weights"""
    # Assume weights are similar when calculating error bars
    assert len(Data) == len(Weights), "Data and Weights size must match!"
    assert len(Weights) > 0, "Data is empty!"
    Z = np.sum(Weights)
    Avg = sum(Data)/sum(Weights)
    if len(Data) > 1:
        Var = sum((d/norm - Avg) ** 2*norm/Z for (d, norm)
                  in zip(Data, Weights))
        return Avg, np.sqrt(Var/(len(Data)-1))
    else:
        return Avg, Var*0.0


def LoadFile(Folder, FileName, Shape):
    Step = []
    Norm = []
    Data = []

    for f in getListOfFiles(Folder):
        if re.search(FileName, f):
            # print "Loading ", f
            try:
                with open(f, "r") as file:
                    Step.append(int(file.readline().split(":")[1]))
                    Norm.append(float(file.readline().split(":")[1]))

                assert len(Norm) == len(Data) + \
                    1, "size of Data and Norm must be the same!"

                Data.append(np.loadtxt(f).reshape(Shape))
            except Exception as e:
                print "Failed to load {0}".format(f)
                print str(e)

    return Data, Norm, Step


def ErrorPlot(p, x, d, color='k', marker='s', label=None, size=4, shift=False):
    p.plot(x, d, marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)


ColorList = ['k', 'r', 'b', 'g', 'm', 'c', 'navy',
             'y', 'cyan', 'darkgreen', 'violet', 'lime', 'purple']
ColorList = ColorList*40


if __name__ == '__main__':
    Para = param()

    dirName = "./Beta40_rs4.0_lambda0.2301"

    for elem in getListOfFiles(dirName):
        print(elem)
