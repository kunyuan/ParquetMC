from numpy import *
import math


def Polarisi(q, w, EF):
    """ Polarization P(q,iW) on imaginary axis. Note that the result is real.
        It works on arrays of frequency, i.e., w can be array of bosonic Matsubara points.
    """
    kF = sqrt(EF)
    q2 = q**2
    kFq = 2*kF*q
    D = 1./(8.*kF*q)

    if type(w) == ndarray:
        res = zeros(len(w), dtype=float)

        # careful for small q or large w
        # for w[i>=iw_start] we should use power expansion
        is_w_large = w > 20*(q2+kFq)
        # if this was newer true, is_w_large contains only False => iw_start=len(w)
        iw_start = len(w)
        # If at least the last frequency is larger than the cutoff, we can find the index
        if is_w_large[-1]:
            iw_start = argmax(is_w_large)

        # if w < cutoff use exact expression
        iw = w[:iw_start]*1j
        wmq2 = iw-q2
        wpq2 = iw+q2

        C1 = log(wmq2-kFq)-log(wmq2+kFq)
        C2 = log(wpq2-kFq)-log(wpq2+kFq)
        res[:iw_start] = real(-kF/(4*pi**2) * (1. - D *
                                               (wmq2**2/q**2-4*EF)*C1 + D*(wpq2**2/q**2-4*EF)*C2))
        # if w < cutoff use proper power expansion
        b2 = q2 * (q2 + 12./5. * EF)  # b2==b^2
        c = 2*EF*kF*q2/(3*pi**2)
        res[iw_start:] = -c/(w[iw_start:]**2 + b2)
    else:
        # careful for small q or large w
        if w == 0.0:
            Nf = kF/2.0/pi/pi
            x = q/2.0/kF
            if abs(x-1.0) < 1.0e-8:
                return -0.5*Nf
            else:
                return -(0.5-(x*x-1)/4.0/x*log(abs((1+x)/(1-x))))*Nf

        if w <= 20*(q2+kFq):
            iw = w*1j
            wmq2 = iw-q2
            wpq2 = iw+q2
            C1 = log(wmq2-kFq)-log(wmq2+kFq)
            C2 = log(wpq2-kFq)-log(wpq2+kFq)
            res = real(-kF/(4*pi**2) * (1. - D*(wmq2**2/q**2-4*EF)
                                        * C1 + D*(wpq2**2/q**2-4*EF)*C2))
        else:
            b2 = q2 * (q2 + 12./5. * EF)
            c = 2*EF*kF*q2/(3*pi**2)
            res = -c/(w**2 + b2)
    return res
