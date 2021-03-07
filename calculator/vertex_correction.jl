using QuantumStatistics
using DelimitedFiles
using Cuba
using NPZ

mass2 = 1
rs = 1
kF = (9π / 4)^(1 / 3) / rs
EF = kF^2
Nf = kF / 2 / π^2
β = 25 / EF
Nt = 64
Nk = 64

T = Grid.tau(β, 1.0 / EF, Nt)
Q = Grid.boseK(kF, 6kF, kF / 8.0, Nk)

open("tgrid.dat", "w") do io
    writedlm(io, T.grid)
end
open("kgrid.dat", "w") do io
    writedlm(io, Q.grid)
end

Rs = npzread("Rs.npy")
Ra = npzread("Ra.npy")

# println(Rs[1, :])

# println(T.grid)
# println(Q.grid)

function Interaction(q, t1, t2)
    t = abs(t1 - t2)
    
    return 8π / (mass2) / β
end

function integrand1(x, f)
    k = x[1] * kF
    t1, t2, t3 = x[2], x[3], x[4]
    t0 = 0.0
    ϵ = (k^2 - kF^2) * β
    g1 = Spectral.kernelFermiT(t3 - t1, ϵ)
    g2 = Spectral.kernelFermiT(t1 - t2, ϵ)
    w1 = Interaction(k, t2, t3)
    w2 = Interaction(0.0, t0, t1)
    phase = exp((t3 - t2) * π) / (2π)^3
    f[1] = g1 * g2 * w1 * w2 * 4π * k^2 * kF * phase * β^3
    # return * 4π * k^2 * kF * phase * β^2
end

function integrand2(x, f)
    k = kF + x[1] / (1 - x[1])
    t1, t2, t3 = x[2], x[3], x[4]
    # t1, t2 = x[2], x[3], x[4]
    t0 = 0.0
    ϵ = (k^2 - kF^2) * β
    # g1 = Spectral.kernelFermiT(t2 - t0, ϵ)
    # g2 = Spectral.kernelFermiT(t0 - t1, ϵ)
    # w = Interaction(k, t1, t2)
    # phase = exp((t2 - t1) * π) / (2π)^3
    # f[1] = g1 * g2 * w * 4π * k^2 * phase / (1 - x[1])^2 * β^3
    g1 = Spectral.kernelFermiT(t3 - t1, ϵ)
    g2 = Spectral.kernelFermiT(t1 - t2, ϵ)
    w1 = Interaction(k, t2, t3)
    w2 = Interaction(0.0, t0, t1)
    phase = exp((t3 - t2) * π) / (2π)^3
    f[1] = g1 * g2 * w1 * w2 * 4π * k^2 * phase / (1 - x[1])^2 * β^3
    # f[1] = g1 * g2 * w1 * w2 * 4π * k^2 * kF * phase * β^3
end

result1, err1 = Cuba.cuhre(integrand1, 4, 1, atol=1e-12, rtol=1e-10);
result2, err2 = Cuba.cuhre(integrand2, 4, 1, atol=1e-12, rtol=1e-10);

println(" Result of Cuba: ", result1[1], " ± ", err1[1])
println(" Result of Cuba: ", result2[1], " ± ", err2[1])

result1, err1 = Cuba.vegas(integrand1, 4, 1, atol=1e-12, rtol=1e-10);
result2, err2 = Cuba.vegas(integrand2, 4, 1, atol=1e-12, rtol=1e-10);

println(" Result of Cuba: ", result1[1], " ± ", err1[1])
println(" Result of Cuba: ", result2[1], " ± ", err2[1])


