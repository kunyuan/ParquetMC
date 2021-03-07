using PyCall
using QuantumStatistics
# math = pyimport("math")
# a = math.sin(math.pi / 4)
# println(a)

pushfirst!(PyVector(pyimport("sys")."path"), "") # add the current directory into the package serach path
# println(PyCall.libpython)

IO = pyimport("utility.IO")
KO = pyimport("KOinteraction")

Para = IO.param("./")
println(Para)

T = Grid.tau(Para.Beta, 1.0 / Para.EF, Para.TauGridSize)
Q = Grid.boseK(Para.kF, Para.MaxExtMom, Para.kF / 8.0, Para.MomGridSize)

println(T.grid)
println(Q.grid)

KO.InterTau(T.grid, Q.grid, Para)

# plt = pyimport("matplotlib.pyplot")
# x = range(0;stop=2 * pi,length=1000); y = sin.(3 * x + 4 * cos.(2 * x));
# plt.plot(x, y, color="red", linewidth=2.0, linestyle="--")
# plt.show()

# ########### Plot Polarization in Tau ################
# plt.figure()
# for qi, q in enumerate(Kbose.grid[:10]):
#     Errorbar(T.grid, dRsT[qi, :], label=f"{q}")
# plt.title("dRs")
# plt.legend()
# plt.show()

# plt.figure()
# for qi, q in enumerate(Kbose.grid[:10]):
#     Errorbar(T.grid, dRaT[qi, :], label=f"{q}")
# plt.title("dRa")
# plt.legend()
# plt.show()


