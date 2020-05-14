import numpy as np
from utility import *

TauGridSize = 1024
MomGridSize = 64
AngGridSize = 64
GridFile = "grid.data"


def BuildTauGrid(Para, GridSize):
    # non-unifor grid
    # Size = GridSize-1
    # arr = np.array(range(Size+1))
    # TauGrid = arr*1.0
    # c = Para.Beta*Para.EF/3.0
    # for i in range(Size/2):
    #     TauGrid[i] = Para.Beta/2.0 * \
    #         np.exp(c*(float(i)/Size-0.5)) - \
    #         Para.Beta/2.0*np.exp(-c*0.5)+1.0e-8
    # for i in range(Size/2, Size+1):
    #     TauGrid[i] = Para.Beta-Para.Beta/2.0 * \
    #         np.exp(c*(0.5-float(i)/Size))+Para.Beta/2.0 * \
    #         np.exp(-c*0.5)-1.0e-8
    # uniform grid
    arr = np.array(range(GridSize))
    TauGrid = arr*Para.Beta/(GridSize-1)
    TauGrid[0] += 1.0e-8
    TauGrid[-1] -= 1.0e-8
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

    TauGrid = BuildTauGrid(Para, TauGridSize)
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
