#!/usr/bin/env python3
import grid
K = grid.FermiK()
K.build(1.0, 3.0, 8, 1.0)
print(K.str())
T = grid.Tau()
T.build(10.0, 16, 3.0)
print(T.str())
Kbose = grid.BoseK()
Kbose.build(1.0, 3.0, 8, 1.0)
print(Kbose.str())
angle = grid.Uniform()
angle.build([0.0, 1.0], 8)
print(angle.str())
