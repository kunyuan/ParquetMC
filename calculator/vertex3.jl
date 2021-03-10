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

results = zeros((3, N))
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

for (qi, q0) in enumerate(LinRange(0.0, 3 * kF, N))
    println("calculating $q0:")

    function vertex3(k, t)
    # bosonic in: (0, 0, q0), fermionic in: (0, 0, kF-q/2), fermionic out: (0, 0, kF+q/2)
        k1 = k # first Green's function
        k2 = @SVector [k[1], k[2], k[3] + q0] # second Green's function
        q = @SVector [k[1], k[2], k[3] - (kF - q0 / 2)]

        ϵ1, ϵ2 = (dot(k1, k1) - kF^2) * β, (dot(k2, k2) - kF^2) * β
        t0 = 0.0
        g1 = Spectral.kernelFermiT(t0 - t[1], ϵ1)
        g0 = Spectral.kernelFermiT(t[1] - t0, ϵ2)
        g2 = Spectral.kernelFermiT(t[2] - t0, ϵ2)

        qd = sqrt(dot(q, q))
        dt = abs(t[2] - t[1]) * β

    # direct part of retared interaction, (t1, t1, t2, t2)
        if (qd <= 1.0e-4)
            wd = 0.0
        else
            wd = -Interpolate.linear2D(Rs, Q, T, qd, dt) * 8π / (qd^2 + mass2)
        end

    # bare interaction part, equal time (t1, t1, t1, t1)
        vd = -8π / (qd^2 + mass2 + 8π * Nf * lindhard(qd / 2.0 / kF)) / β
        vd -= wd

        phase0 = cos(π * (t[1] + t[1])) / (2π)^3
        phase1 = cos(π * (t[1] + t[2])) / (2π)^3

        gamma = g1 * g0 * vd * phase0 + g1 * g2 * wd * phase1

        return gamma
    end

    function counterterm(k, t)
        ϵ1, ϵ2 = (dot(k, k) - kF^2) * β, (dot(k, k) - kF^2) * β
        t0 = 0.0
        g0 = Spectral.kernelFermiT(t[1] - t0, ϵ2)
        g1 = Spectral.kernelFermiT(t0 - t[1], ϵ1)
        g2 = Spectral.kernelFermiT(t[2] - t0, ϵ2)
        q = @SVector [k[1], k[2], k[3] + kF]

        kr = k / sqrt(dot(k, k)) * kF

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

    result1, err1 = Cuba.cuhre(integrand1, 5, 1, atol=1e-12, rtol=1e-10);
    result2, err2 = Cuba.cuhre(integrand2, 5, 1, atol=1e-12, rtol=1e-10);

# println(" Dir: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
    println(" q0=$q0 (cuhre): ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
    results[1, qi] = q0
    results[2, qi] = result1[1] + result2[1]

    result1, err1 = Cuba.vegas(integrand1, 5, 1, atol=1e-12, rtol=1e-10);
    result2, err2 = Cuba.vegas(integrand2, 5, 1, atol=1e-12, rtol=1e-10);

# println(" Result of Cuba: ", result1[1], " ± ", err1[1])
# println(" Result of Cuba: ", result2[1], " ± ", err2[1])
# println(" Total: ", result1[1] + result2[1], " ± ", err1[1] + err2[1])
    println(" q0=$q0 (vegas): ", result1[1] + result2[1], " ± ", err1[1] + err2[1])

    results[3, qi] = result1[1] + result2[1]
end

for qi in 1:N
    println("$(results[1, qi] / kF)    $(results[2, qi])    $(results[3, qi])")
end
