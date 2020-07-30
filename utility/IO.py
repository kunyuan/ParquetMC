from utility.color import *
import numpy as np
import glob
import re
import sys
import os


def GetLine(file):
    while True:
        line = file.readline().strip()
        if len(line) > 0 and line[0] != "#":
            return line


def getListOfFiles(dirName):
    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(dirName):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    return listOfFiles


class param:
    # Order, Beta, Rs, Mass2, Lambda, Charge2, TotalStep = [None, ]*7
    # kF, Nf, EF, Bubble = [0.0, ]*4
    def __init__(self):
        self.DataFolder = "Data"
        self.InputFile = "parameter"

        with open(self.InputFile, "r") as file:
            file.readline()  # comment line
            para = GetLine(file).split(",")
            self.Order = int(para[0])
            self.Beta = float(para[1])
            self.Rs = float(para[2])
            self.Mass2 = float(para[3])
            self.Lambda = float(para[4])
            self.Charge2 = float(para[5])
            self.Dim = int(para[6])
            self.Spin = int(para[7])
            self.TotalStep = int(para[8])

            grid = GetLine(file).split(",")
            self.TauGridSize = int(grid[0])
            self.MomGridSize = int(grid[1])
            self.AngGridSize = int(grid[2])
            self.MaxExtMom = float(grid[3])

            timer = GetLine(file).split(",")
            self.PrintTimer = int(timer[0])
            self.SaveTimer = int(timer[1])
            self.ReWeightTimer = int(timer[2])
            self.MessageTimer = int(timer[3])
            self.CollectionTimer = int(timer[4])

        if self.Dim == 3:
            self.kF = (9.0*np.pi/4.0)**(1.0/3.0)/self.Rs
            self.Nf = self.kF/4.0/np.pi**2*self.Spin
        elif self.Dim == 2:
            self.kF = np.sqrt(2.0)/self.Rs  # 2D
            self.Nf = 1.0/4.0/np.pi*self.Spin
        else:
            print(f"Not Implemented for Dimension {self.Dim}")
            sys.exit(0)

        self.EF = self.kF**2
        self.Beta /= self.EF
        self.MaxExtMom *= self.kF

        print(yellow("Parameters:"))
        print(f"Rs={self.Rs}, kF={self.kF}, EF={self.EF}, Beta={self.Beta}, Mass2={self.Mass2}, Lambda={self.Lambda}, Dim={self.Dim}, Spin={self.Spin}\n")
        print(yellow("Grid Information:"))
        print(
            f"TauSize={self.TauGridSize}, MomSize={self.MomGridSize}, AngleSize={self.AngGridSize}, MaxExtMom={self.MaxExtMom}")

        print(yellow("Timer Information:"))
        print(f"Print={self.PrintTimer}, Save={self.SaveTimer}, ReWeight={self.ReWeightTimer}, Message={self.MessageTimer}, Collection={self.CollectionTimer}")

# def Estimate(Data, Weights, axis=0):
#     """ Return Mean and Error  with given weights"""
#     # Assume weights are similar when calculating error bars
#     Weights = np.array(Weights)
#     Num = len(Weights)
#     assert Num > 0, "Data is empty!"
#     assert Data.shape[0] == Num, "Data and Weights size must match!"
#     Avg = np.average(Data, weights=Weights, axis=0)
#     Var = np.average((Data-Avg)**2, weights=Weights, axis=0)
#     Err = np.sqrt(Var/(Num-1)) if Num > 1 else 0.0
#     return Avg, Err


def Estimate(Data, Weights, operation=None):
    """ Return Mean and Error  with given weights"""
    # Assume weights are similar when calculating error bars
    assert len(Data) == len(Weights), "Data and Weights size must match!"
    assert len(Weights) > 0, "Data is empty!"

    if operation is not None:
        Data = [operation(d) for d in Data]

    Z = np.sum(Weights)
    Avg = sum(Data)/sum(Weights)
    if len(Data) > 1:
        Var = sum((d/norm - Avg) ** 2*norm/Z for (d, norm)
                  in zip(Data, Weights))
        return Avg, np.sqrt(Var/(len(Data)-1))
    else:
        return Avg, Var*0.0


def LoadFile(Folder, FileName, shape=None):
    Step, Norm, Data = [], [], []
    Grid = {}

    for f in getListOfFiles(Folder):
        if re.search(FileName, f):
            print("Loading ", f)
            try:
                with open(f, "r") as file:
                    Step.append(int(file.readline().split(":")[1]))
                    Norm.append(float(file.readline().split(":")[1]))
                    while True:
                        g = file.readline().split(":")
                        if g[0].find("Grid") != -1:
                            key = g[0].strip(" #")
                            Grid[key] = np.fromstring(g[1], sep=' ')
                        else:
                            break

                assert len(Norm) == len(Data) + \
                    1, "size of Data and Norm must be the same!"

                if shape == None:
                    Data.append(np.loadtxt(f, dtype=float))
                else:
                    Data.append(np.loadtxt(f).reshape(shape))

            except Exception as e:
                print(f"Failed to load {f}")
                print(str(e))

    return Data, Norm, Step, Grid


if __name__ == '__main__':
    Para = param()

    dirName = "../Data"
    filename = "sigma_pid[0-9]+.dat"

    LoadFile(dirName, filename)
    # for elem in getListOfFiles(dirName):
    #     print(elem)
