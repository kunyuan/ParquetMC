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

IsF = false
INL, OUTL, INR, OUTR = 1, 2, 3, 4
DI, EX = 1, 2

k0 = @SVector [0.0, 0.0, kF]

legK = [k0, k0, k0, k0]

Rs = npzread("Rs.npy") # dimension: k, t
Ra = npzread("Ra.npy") # dimension: k, t

# println(Rs[1, :])

# println(T.grid)
# println(Q.grid)
function lindhard(x)
    if (abs(x) < 1.0e-4)
        return 1.0
    elseif (abs(x - 1.0) < 1.0e-4)
        return 0.5
    else
        return 0.5 - (x^2 - 1) / 4.0 / x * log(abs((1 + x) / (1 - x)))
    end
end

function phase(tInL, tOutL, tInR, tOutR)
    # return 1.0;
    if (IsF)
        return cos(π * ((tInL + tOutL) - (tInR + tOutR)));
    else
        return cos(π * ((tInL - tOutL) + (tInR - tOutR)))
    end
end

function vertex3(k, t)
    ϵ1, ϵ2 = (dot(k, k) - kF^2) * β, (dot(k, k) - kF^2) * β
    t0 = 0.0
    g0 = Spectral.kernelFermiT(t[1] - t0, ϵ2)
    g1 = Spectral.kernelFermiT(t0 - t[1], ϵ1)
    g2 = Spectral.kernelFermiT(t[2] - t0, ϵ2)
    q = @SVector [k[1], k[2], k[3] + kF]

    qd = sqrt(dot(q, q))
    dt = abs(t[2] - t[1]) * β

    # bare interaction part, equal time (t1, t1, t1, t1)

    # direct part of retared interaction, (t1, t1, t2, t2)
    if (qd <= 1.0e-4)
        wd = 0.0
    else
        wd = -Interpolate.linear2D(Rs, Q, T, qd, dt) * 8π / (qd^2 + mass2)
    end

    vd = -8π / (qd^2 + mass2 + 8π * Nf * lindhard(qd / 2.0 / kF)) / β
    vd -= wd

    # vd = -8π / (qd^2 + mass2) / β

    phase0 = cos(π * (t[1] - t[1])) / (2π)^3
    phase1 = cos(π * (t[1] - t[2])) / (2π)^3

    # phase0 = cos(π * (t[1] - t[1])) / (2π)^3
    # phase1 = cos(π * (t[1] - t[2])) / (2π)^3

    # gamma = (g1 * g2 + g1 * g0) * phase

    gamma = g1 * g0 * vd * phase0 + g1 * g2 * wd * phase1
    # gamma = g1 * g0 * vd * phase0

    return gamma
end

function integrand1(x, f)
    k, θ, ϕ = x[1] * kF, x[2] * π, x[3] * 2π
    k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
    # t = @SVector [0.0, x[4], x[5], x[6]]
    t = @SVector [x[4], x[5]]
    factor = kF * 2π^2 * β^2 * k^2 * sin(θ) 
    f[1] = vertex3(k1, t) * factor * Nf 
    # return * 4π * k^2 * kF * phase * β^2
end

function integrand2(x, f)
    k, θ, ϕ = kF + x[1] / (1 - x[1]), x[2] * π, x[3] * 2π
    k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
    # t = @SVector [0.0, x[4], x[5], x[6]]
    t = @SVector [x[4], x[5]]
    factor = 2π^2 * β^2 * k^2 * sin(θ) / (1 - x[1])^2
    f[1] = vertex3(k1, t) * factor * Nf 
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
#     w1 = 8π / (q^2 + mass2) / β
#     # w1 = bare(q, t1, t2)
#     # w2 = Interaction(0.0, t0, t1)
#     phase = cos((t2 - t2) * π) / (2π)^3 
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
#     w1 = 8π / (q^2 + mass2) / β
#     # w2 = Interaction(0.0, t0, t1)
#     phase = cos((t2 - t2) * π) / (2π)^3 
#     factor = 2π^2 * β^2 / (1 - x[1])^2 
#     f[1] = w1 * g1 * g2 * k^2 * sin(θ) * phase  * factor
#     # f[1] = g1 * g2 * w1 * w2 * 4π * k^2 * kF * phase * β^3

# end

result1, err1 = Cuba.cuhre(integrand1, 5, 1, atol=1e-12, rtol=1e-10);
result2, err2 = Cuba.cuhre(integrand2, 5, 1, atol=1e-12, rtol=1e-10);

# println(" Dir: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
println(" Dir: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])

result1, err1 = Cuba.vegas(integrand1, 5, 1, atol=1e-12, rtol=1e-10);

# println(" Result of Cuba: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
# println(" Total: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println(" Dir: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])


