import numpy as np
from utility import *

TauGridSize = 128
MomGridSize = 64
AngGridSize = 64
GridFile = "grid.data"


def BuildTauGrid(Beta, GridSize):
    arr = np.array(range(GridSize))
    TauGrid = arr*Beta/GridSize+1.0e-8
    return TauGrid


def BuildMomGrid(MaxK, GridSize):
    arr = np.array(range(GridSize))
    KGrid = arr*MaxK/GridSize+1.0e-8
    return KGrid


def BuildAngleGrid(GridSize):
    arr = np.array(range(GridSize))
    AngGrid = arr*2.0/GridSize-1.0
    return AngGrid


if __name__ == "__main__":
    Para = param()

    TauGrid = BuildTauGrid(Para.Beta, TauGridSize)
    MomGrid = BuildMomGrid(Para.MaxExtMom, MomGridSize)
    AngGrid = BuildAngleGrid(AngGridSize)

    with open(GridFile, "w") as f:
        f.writelines("{0} #TauGrid\n".format(TauGridSize))
        for t in TauGrid:
            f.write("{0} ".format(t))
        f.write("\n\n")

        f.writelines("{0} #MomGrid\n".format(MomGridSize))
        for k in MomGrid:
            f.write("{0} ".format(k))
        f.write("\n\n")

        f.writelines("{0} #AngGrid\n".format(AngGridSize))
        for a in AngGrid:
            f.write("{0} ".format(a))
        f.write("\n\n")

    pass
