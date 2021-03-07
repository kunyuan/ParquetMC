using QuantumStatistics
using DelimitedFiles

include("parameter.jl")

T = Grid.tau(β, 1.0 / EF, Nt)
Q = Grid.boseK(kF, 6kF, kF / 8.0, Nk)

open("tgrid.dat", "w") do io
    writedlm(io, T.grid)
end
open("kgrid.dat", "w") do io
    writedlm(io, Q.grid)
end
# println("tgrid: ", T.grid)
# println("kgrid: ", Q.grid)