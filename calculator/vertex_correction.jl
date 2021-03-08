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

function interaction(qd, qe, t1, t2)
    qd = sqrt(dot(qd, qd))
    qe = sqrt(dot(qe, qe))
    dt = abs(t2 - t1) * β

    # bare interaction part, equal time (t1, t1, t1, t1)
    vd = -8π / (qd^2 + mass2) / β
    ve = 8π / (qe^2 + mass2) / β

    # direct part of retared interaction, (t1, t1, t2, t2)
    if (qd <= Q.grid[1])
        wd = -Interpolate.linear2D(Rs, Q, T, 1.0e-6, dt)
    else
        wd = -Interpolate.linear2D(Rs, Q, T, qd, dt)
    end

    # exchange part of retared interaction (t1, t2, t2, t1)
    if (qe <= Q.grid[1])
        we = Interpolate.linear2D(Rs, Q, T, 1.0e-6, dt)
    else
        we = Interpolate.linear2D(Rs, Q, T, qe, dt)
    end

    # vd, ve = 0.0, 0.0
    # wd = 1.0 / (qd * qd + 1)
    # we = 1.0 / (qe * qe + 1)

    # println(wd, ", ", vd)

    # ve, we = 0.0, 0.0

    return vd, ve, wd, we
    # return 1 / β
end

function phase(tInL, tOutL, tInR, tOutR)
    # return 1.0;
    if (IsF)
        return cos(π * ((tInL + tOutL) - (tInR + tOutR)));
    else
        return cos(π * ((tInL - tOutL) + (tInR - tOutR)))
    end
end

function Tchannel(legK, t, k1)
    qd = legK[INL] - legK[OUTL]
    k2 = k1 - qd
    vld, vle, wld, wle = interaction(qd, legK[INL] - k1, t[1], t[2])
    vrd, vre, wrd, wre = interaction(qd, legK[INR] - k2, t[3], t[4])

    # println("$vld, $vrd, $wld, $wrd")

    ϵ1, ϵ2 = (dot(k1, k1) - kF^2) * β, (dot(k2, k2) - kF^2) * β
    wd, we = 0.0, 0.0

    gu13 = Spectral.kernelFermiT(t[3] - t[1], ϵ1)
    gu23 = Spectral.kernelFermiT(t[3] - t[2], ϵ1)

    gd31 = Spectral.kernelFermiT(t[1] - t[3], ϵ2)
    gd41 = Spectral.kernelFermiT(t[1] - t[4], ϵ2)
    gd32 = Spectral.kernelFermiT(t[2] - t[3], ϵ2)
    gd42 = Spectral.kernelFermiT(t[2] - t[4], ϵ2)

    # # v(1111)*G(k1, 13)*G(k2, 31)*v(3333) ==> (1133)
    G = gu13 * gd31 / (2π)^3 * phase(t[1], t[1], t[3], t[3])
    wd += G * (SPIN * vld * vrd + vld * vre + vle * vrd)
    we += G * (vle * vre)

    # v(1111)*G(k1, 13)*G(k2, 31)*Wd(3344) ==> (1144)
    G = gu13 * gd31 / (2π)^3 * phase(t[1], t[1], t[4], t[4])
    wd += G * (SPIN * vld * wrd + vle * wrd)
    we += 0.0

    # Wd(1122)*G(k1, 23)*G(k2, 32)*v(3333) ==> (1133)
    G = gu23 * gd32 / (2π)^3 * phase(t[1], t[1], t[3], t[3])
    wd += G * (SPIN * wld * vrd + wld * vre)
    we += 0.0

    # Wd(1122)*G(k1, 23)*G(k2, 32)*Wd(3344) ==> (1144)
    G = gu23 * gd32 / (2π)^3 * phase(t[1], t[1], t[4], t[4])
    wd += G * (SPIN * wld * wrd + vle * wrd)
    we += 0.0

    # v(1111)*G(k1, 13)*G(k2, 41)*We(3443) ==> (1143)
    G = gu13 * gd41 / (2π)^3 * phase(t[1], t[1], t[4], t[3])
    wd += G * (vld * wre)
    we += G * (vle * wre)

    # Wd(1122)*G(k1, 24)*G(k2, 32)*We(3443) ==> (1143)
    G = gu23 * gd42 / (2π)^3 * phase(t[1], t[1], t[4], t[3])
    wd += G * (wld * wre)
    we += 0.0

    # # We(1221)*G(k1, 13)*G(k2, 32)*v(3333) ==> (1233)
    G = gu13 * gd32 / (2π)^3 * phase(t[1], t[2], t[3], t[3])
    wd += G * (wle * vrd)
    we += G * (wle * vre)

    # We(1221)*G(k1, 13)*G(k2, 32)*Wd(3344) ==> (1244)
    G = gu13 * gd32 / (2π)^3 * phase(t[1], t[2], t[4], t[4])
    wd += G * (wle * wrd)
    we += 0.0

    # We(1221)*G(k1, 13)*G(k2, 42)*We(3443) ==> (1243)
    G = gu13 * gd42 / (2π)^3 * phase(t[1], t[2], t[4], t[3])
    wd += 0.0
    we += G * (wle * wre)

    return @SVector [wd, we]
end

function integrand1(x, f)
    k, θ, ϕ = x[1] * kF, x[2] * π, x[3] * 2π
    k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
    # t = @SVector [0.0, x[4], x[5], x[6]]
    t = @SVector [x[4], x[5], x[6], x[7]]
    factor = kF * 2π^2 * β^3 * k^2 * sin(θ) 
    f[1], f[2] = Tchannel(legK, t, k1) * factor * Nf 
    # return * 4π * k^2 * kF * phase * β^2
end

function integrand2(x, f)
    k, θ, ϕ = kF + x[1] / (1 - x[1]), x[2] * π, x[3] * 2π
    k1 = @SVector [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
    # t = @SVector [0.0, x[4], x[5], x[6]]
    t = @SVector [x[4], x[5], x[6], x[7]]
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

result1, err1 = Cuba.cuhre(integrand1, 7, 2, atol=1e-12, rtol=1e-10);
result2, err2 = Cuba.cuhre(integrand2, 7, 2, atol=1e-12, rtol=1e-10);

# println(" Dir: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
println(" Dir: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println(" Ex : ", result1[2] + result2[2], " ± ", err1[2] + err2[2])

result1, err1 = Cuba.vegas(integrand1, 7, 2, atol=1e-12, rtol=1e-10);
result2, err2 = Cuba.vegas(integrand2, 7, 2, atol=1e-12, rtol=1e-10);

# println(" Result of Cuba: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
# println(" Total: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println(" Dir: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
println(" Ex : ", result1[2] + result2[2], " ± ", err1[2] + err2[2])


