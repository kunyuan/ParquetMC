#!/usr/bin/env python3
import grid
K = grid.FermiK()
K.build(1.0, 3.0, 8, 1.0)  # double kF, double maxK, int _size, double scale
print(K.str())
T = grid.Tau()
T.build(10.0, 16, 3.0)  # double beta, int _size, double scale
print(T.str())
Kbose = grid.BoseK()
# double kF, double maxK, int _size, double scale
Kbose.build(1.0, 3.0, 8, 1.0)
print(Kbose.str())
angle = grid.Uniform()
angle.build([0.0, 1.0], 8)  # std::array<double, 2> bounds, int _size
print(angle.str())
