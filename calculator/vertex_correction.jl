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

INL, OUTL, INR, OUTR = 1, 2, 3, 4

k0 = @SVector [0.0, 0.0, kF]

legK = [k0, k0, k0, k0]

Rs = npzread("Rs.npy")
Ra = npzread("Ra.npy")

# println(Rs[1, :])

# println(T.grid)
# println(Q.grid)

function bare(qd, qe, t1, t2)
    qd2 = dot(qd, qd)
    qe2 = dot(qe, qe)
    wd = -8π / (qd2 + mass2) / β
    we = 8π / (qe2 + mass2) / β
    return wd, we
    # return 1 / β
end

function Tchannel(legK, t, k1)
    qd = legK[INL] - legK[OUTL]
    k2 = qd - k1
    g1 = Spectral.kernelFermiT(t[3] - t[1], (dot(k1, k1) - kF^2) * β)
    g2 = Spectral.kernelFermiT(t[1] - t[3], (dot(k2, k2) - kF^2) * β)
    wld, wle = bare(qd, legK[INL] - k1, t[1], t[2])
    wrd, wre = bare(qd, legK[INR] - k2, t[3], t[4])

    # vd = SPIN * wld * wrd + wld * wre + wle * wrd
    vd = SPIN * wld * wrd
    ve = wle * wre

    G = g1 * g2
    phase = exp((t[3] - t[3]) * π) / (2π)^3 * (-1)

    wd, we = phase * G * vd, phase * G * ve
    return @SVector [wd, we]
end

function integrand1(x, f)
    k, θ, ϕ = x[1] * kF, x[2] * π, x[3] * 2π
    k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
    t = @SVector [0.0, x[4], x[5], x[6]]
    factor = kF * 2π^2 * β^3 * k^2 * sin(θ) 
    f[1], f[2] = Tchannel(legK, t, k1) * factor * Nf
    # return * 4π * k^2 * kF * phase * β^2
end

function integrand2(x, f)
    k, θ, ϕ = kF + x[1] / (1 - x[1]), x[2] * π, x[3] * 2π
    k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
    t = @SVector [0.0, x[4], x[5], x[6]]
    factor = 2π^2 * β^3 * k^2 * sin(θ) / (1 - x[1])^2
    f[1], f[2] = Tchannel(legK, t, k1) * factor * Nf
    # return * 4π * k^2 * kF * phase * β^2
end


# function integrand1(x, f)
#     k, θ, ϕ = x[1] * kF, x[2] * π, x[3] * 2π
#     kx, ky, kz = k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)
#     t0, t1, t2 = 0.0, x[4], x[5]

#     ϵ = (k^2 - kF^2) * β
#     g1 = Spectral.kernelFermiT(t1 - t0, ϵ)
#     g2 = Spectral.kernelFermiT(t0 - t1, ϵ)
#     q = sqrt((kz + kF)^2 + kx^2 + ky^2)
#     w1 = bare(q, t1, t2)
#     # w2 = Interaction(0.0, t0, t1)
#     phase = exp((t2 - t2) * π) / (2π)^3 
#     factor = kF * 2π^2 * β^2 
#     f[1] = w1 * g1 * g2 * k^2 * sin(θ) *  phase * factor

#     # return * 4π * k^2 * kF * phase * β^2
# end

# function integrand2(x, f)
#     k, θ, ϕ = kF + x[1] / (1 - x[1]), x[2] * π, x[3] * 2π
#     kx, ky, kz = k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)
#     t0, t1, t2 = 0.0, x[4], x[5]

#     ϵ = (k^2 - kF^2) * β
#     g1 = Spectral.kernelFermiT(t1 - t0, ϵ)
#     g2 = Spectral.kernelFermiT(t0 - t1, ϵ)
#     q = sqrt((kz + kF)^2 + kx^2 + ky^2)
#     w1 = bare(q, t1, t2)
#     # w2 = Interaction(0.0, t0, t1)
#     phase = exp((t2 - t2) * π) / (2π)^3 
#     factor = 2π^2 * β^2 / (1 - x[1])^2 
#     f[1] = w1 * g1 * g2 * k^2 * sin(θ) * phase  * factor
#     # f[1] = g1 * g2 * w1 * w2 * 4π * k^2 * kF * phase * β^3

# end

result1, err1 = Cuba.cuhre(integrand1, 6, 2, atol=1e-12, rtol=1e-10);
result2, err2 = Cuba.cuhre(integrand2, 6, 2, atol=1e-12, rtol=1e-10);

# println(" Dir: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
println(" Dir: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println(" Ex : ", result1[2] + result2[2], " ± ", err1[2] + err2[2])

result1, err1 = Cuba.vegas(integrand1, 6, 2, atol=1e-12, rtol=1e-10);
result2, err2 = Cuba.vegas(integrand2, 6, 2, atol=1e-12, rtol=1e-10);

# println(" Result of Cuba: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
# println(" Total: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println(" Dir: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println(" Ex : ", result1[2] + result2[2], " ± ", err1[2] + err2[2])


