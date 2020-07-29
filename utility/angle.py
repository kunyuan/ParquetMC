from scipy.special import eval_legendre
from scipy import integrate
import sys
import os
import numpy as np
import unittest


def mult_along_axis(A, B, axis):
    """
    return A[..., i, ...]*B[i],
    where B[i] is broad casted to all elements of A[..., i, ...]
    """

    A = np.array(A)
    B = np.array(B)

    # shape check
    if axis >= A.ndim:
        raise ValueError(r"{axis} is larger than {A.ndim}")
    if A.shape[axis] != B.size:
        raise ValueError(
            "'A' and 'B' must have the same length along the given axis")

    # Swap the given axis with the last one, to get the swapped
    #     'shape' tuple here (swapay only returns a view of the
    #     supplied array, so no unneccesary copy)
    shape = np.swapaxes(A, A.ndim-1, axis).shape

    # Broadcast to an array with the shape as above(no new array created)
    B_brc = np.broadcast_to(B, shape)

    # Swap back the axes (again, this is only a view of the same array)
    B_brc = np.swapaxes(B_brc, A.ndim-1, axis)

    return A * B_brc


def LegendreCoeff(Data, AngleGrid, lList, axis=0):
    """
    Calculate the Legendre Coefficients: \frac{1}{2} \int_{-1}^{1} f(x)*P_l (x) dx
    Data: array with at least one angle axis
    AngleGrid: Grid for cos(theta) in [-1, 1]
    lList: list of angular momentum quantum number
    axis: axis to perform angle integration
    """
    Coeff = {}
    for l in lList:
        data = np.copy(Data)
        legendreP = np.array([eval_legendre(l, x) for x in AngleGrid])
        data = mult_along_axis(data, legendreP, axis)
        Coeff[l] = integrate.simps(data, AngleGrid, axis=axis)/2.0
    return Coeff


def AngleFunc(Coeff, AngleGrid):
    """
    f(x)=\sum_0^∞ (2l+1)*C_l*P_l(x)
    """
    data = np.zeros_like(AngleGrid)
    for l in Coeff.keys():
        legendreP = np.array([eval_legendre(l, x) for x in AngleGrid])
        data += (2.0*l+1)*Coeff[l]*legendreP
    return data


class TestAngle(unittest.TestCase):
    def test_LegendreCoff(self):
        grid = np.linspace(-1.0, 1.0, 128)
        l = 3
        def poly(x): return 0.5*(5.0*x**3-3.0*x)
        d = np.array([poly(x) for x in grid])
        data = np.zeros((2, len(grid)))
        data[0, :] = d
        data[1, :] = d*2.0
        coeff = LegendreCoeff(data, grid, [l, ], axis=1)[l]
        coeff *= 2.0*l+1.0
        # print(coeff)
        self.assertTrue(abs(coeff[0]-1.0) < 1.0e-3)
        self.assertTrue(abs(coeff[1]-2.0) < 1.0e-3)

    def test_backforth(self):
        grid = np.linspace(-1.0, 1.0, 128)
        l = range(6)
        y = grid**3
        coeff = LegendreCoeff(y, grid, l)
        yp = AngleFunc(coeff, grid)
        # print(abs(yp-y))
        # print(np.amax(abs(yp-y)))
        self.assertTrue(np.amax(abs(y-yp)) < 1.0e-3)


if __name__ == '__main__':
    unittest.main()
