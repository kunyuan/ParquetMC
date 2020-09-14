from scipy import integrate
from utility import *
from grid import *

XType = "Tau"
# XType = "Mom"

Para = param()
Grid = grid(Para)
Order = range(0, Para.Order+1)
TauGrid = Grid.TauGrid
MomGrid = Grid.MomGrid

TauBinSize=Para.TauGridSize
ExtMomBinSize=Para.MomGridSize
TauBin=TauGrid
ExtMomBin=MomGrid
with open("./Data/f.dat","w") as file:
    for i in range(Para.TauGridSize):
        file.write("{0} ".format(TauGrid[i]))
        file.write("\t")
    file.write("\n")
    for i in range(TauBinSize):
        for k in range(ExtMomBinSize):
            file.write("{0}\t".format( np.exp(-ExtMomBin[k]**2)*(Para.Beta-2*TauBin[i]) ))