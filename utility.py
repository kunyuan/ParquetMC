import os
import sys
import re
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mat
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

##################### Global variables  #############################
Dim = 3
SpinIndex = 2
InputFile = "inlist"
#####################################################################


class param:
    Order, Beta, Rs, Mass2, Lambda, Charge2, TotalStep = [None, ]*7
    kF, Nf, EF, Bubble = [0.0, ]*4

    def __init__(self):

        with open(InputFile, "r") as file:
            para = file.readline().split(" ")
            self.Order = int(para[0])
            self.Beta = float(para[1])
            self.Rs = float(para[2])
            self.Mass2 = float(para[3])
            self.Lambda = float(para[4])
            self.Charge2 = float(para[5])
            self.MaxExtMom = float(para[6])
            self.TotalStep = int(para[7])

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
        print "Rs={0}, kF={1}, EF={2}, Beta={3}".format(
            self.Rs, self.kF, self.EF, self.Beta)


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


def main():
    dirName = "./Beta40_rs4.0_lambda0.2301"

    for elem in getListOfFiles(dirName):
        print(elem)


if __name__ == '__main__':
    main()
