using QuantumStatistics
using DelimitedFiles
using LinearAlgebra
using Cuba
using NPZ
using StaticArrays

# include("parameter.jl")
include("grid.jl")

println("rs=$rs mass2=$mass2 kF=$kF β=$β Nf=$Nf")
P0, err = Diagram.Polar0(1.0e-4, 0, kF, β)
NfT = Diagram.println("Nf(T)=$(real(P0) * SPIN)  ± $(real(err) * SPIN)")

N = 8
IsF = false
INL, OUTL, INR, OUTR = 1, 2, 3, 4
DI, EX = 1, 2

Rs = npzread("Rs.npy") # dimension: k, t
Ra = npzread("Ra.npy") # dimension: k, t

k0 = 2.0kF
w0 = π

function lindhard(x)
    if (abs(x) < 1.0e-4)
        return 1.0
    elseif (abs(x - 1.0) < 1.0e-4)
        return 0.5
    else
        return 0.5 - (x^2 - 1) / 4.0 / x * log(abs((1 + x) / (1 - x)))
    end
end

function sigma(k, t)
    q = @SVector [-k[1], -k[2], -k[3] + k0]
    ϵ = (dot(k, k) - kF^2) * β
    g0 = Spectral.kernelFermiT(-1.0e-6, ϵ)
    g1 = Spectral.kernelFermiT(t, ϵ)

    qd = sqrt(dot(q, q))
    dt = abs(t) * β

    # direct part of retared interaction, (t1, t1, t2, t2)
    if (qd <= 1.0e-4)
        wd = 0.0
    else
        wd = -Interpolate.linear2D(Rs, Q, T, qd, dt) * 8π / (qd^2 + mass2)
    end

    # bare interaction part, equal time (t1, t1, t1, t1)
    vd = -8π / (qd^2 + mass2 + 8π * Nf * lindhard(qd / 2.0 / kF)) / β
    vd -= wd

    # vd = -8π / (qd^2 + mass2 + 8π * Nf) / β
    # vd = -8π / (qd^2 + mass2) / β
    # wd = 0.0

    phase0 = 1 / (2π)^3
    phase1 = cos(w0 * t) / (2π)^3  # real part
    phase2 = sin(w0 * t) / (2π)^3  # imag part

    sigma1 = g0 * vd * phase0 + g1 * wd * phase1
    sigma2 = g1 * wd * phase2

    return sigma1, sigma2
end

function integrand1(x, f)
    k, θ, ϕ = x[1] * 10kF, x[2] * π, x[3] * 2π
    k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
# t = @SVector [0.0, x[4], x[5], x[6]]
    t = x[4]
    factor = 10kF * 2π^2 * β * k^2 * sin(θ) 
    f[1], f[2] = sigma(k1, t)
    f[1] *= factor
    f[2] *= factor
# return * 4π * k^2 * kF * phase * β^2
end

# function integrand2(x, f)
#     k, θ, ϕ = x[1] / (1 - x[1]), x[2] * π, x[3] * 2π
#     k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
# # t = @SVector [0.0, x[4], x[5], x[6]]
#     t = x[4]
#     factor = 2π^2 * β^2 * k^2 * sin(θ) / (1 - x[1])^2
#     f[1], f[2] = sigma(k1, t)
#     f[1] *= factor
#     f[2] *= factor
# # return * 4π * k^2 * kF * phase * β^2
# end

result1, err1 = Cuba.cuhre(integrand1, 4, 2, atol=1e-12, rtol=1e-10);
# result2, err2 = Cuba.cuhre(integrand2, 4, 2, atol=1e-12, rtol=1e-10);
println("(cuhre): real: ", result1[1], " ± ", err1[1])
println("(cuhre): imag: ", result1[2], " ± ", err1[2])

result1, err1 = Cuba.vegas(integrand1, 4, 2, atol=1e-12, rtol=1e-10);
# result2, err2 = Cuba.vegas(integrand2, 4, 2, atol=1e-12, rtol=1e-10);

# println(" Result of Cuba: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
# println(" Total: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println("(vegas): real: ", result1[1], " ± ", err1[1])
println("(vegas): imag: ", result1[2], " ± ", err1[2])
