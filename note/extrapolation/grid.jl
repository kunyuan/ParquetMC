module Grid
using StaticArrays

struct Coeff
    bound::SVector{2,Float64}
    idx::SVector{2,Float64}
    lambda::Float64
    a::Float64
    b::Float64

    function Coeff(b, i, l, dense2sparse)
        bound = b
        idx = i
        lambda = l

        if dense2sparse == false
            bound[1], bound[2] = bound[2], bound[1]
            idx[1], idx[2] = idx[2], idx[1]
            lambda *= -1
        end
        _l1, _l2 = 1.0, exp(lambda * (idx[2] - idx[1]))
        b = (bound[2] - bound[1]) / (_l2 - _l1)
        a = (bound[1] * _l2 - bound[2] * _l1) / (_l2 - _l1)
        return new(bound, idx, lambda, a, b)
    end
end

function _floor(l::Coeff, x)
    @assert l.bound[0] <= x <= l.bound[1] "out of range!"
    pos = l.idx[1] + 1.0 / l.lambda * log((x - l.a) / l.b)
    return Base.floor(Int, pos)
end

function _grid(l::Coeff, idx)
    return l.a + l.b * exp(l.lambda * (Float64(idx) - l.idx[1]))
end

struct LogGrid
    size::Int
    grid::Array{Float64,1}
    segmentEnds::Array{Int,1} # ends of each segments
    coeff::Array{Coeff,1}

    function LogGrid(coeff::Array{Coeff,1}, segmentEnds::Array{Int,1})
        @assert length(segmentEnds) == length(coeff) "Number of coeff and segments don't match!"
        size = segmentEnds[end]
        grid = zeros(Float64, size)
        # println("length= $(length(grid))")
        start = 1
        for (i, seg) in enumerate(segmentEnds)
            for idx = start:seg
                grid[idx] = _grid(coeff[i], idx)
                # println(idx)
            end
            start = seg
        end
        @assert size == length(grid) "length of the grid doesn't match with the size!"
        return new(size, grid, segmentEnds, coeff)
    end
end

struct UniformGrid
    size::Int
    grid::Array{Float64,1}

    function UniformGrid(head, tail, size)
        @assert size > 1 "Size must be large than 1"
        grid = LinRange(head, tail, Int(size))
        return new(Int(size), grid)
    end
end

# cos(theta) grids
function angleGrid(AngGridSize)
    return Grid.UniformGrid(-1.0, 1.0, AngGridSize)
end

# Tau Grid construction
function tauGrid(Beta, Ef, TauGridSize)
    lambda = Beta / Ef / TauGridSize / 3.0
    c1 = Grid.Coeff([0.0, Beta / 2.0], [1.0, TauGridSize / 2 + 0.5], lambda, true)
    c2 = Grid.Coeff([Beta / 2.0, Beta], [TauGridSize / 2 - 0.5, TauGridSize], lambda, false)
    tau = Grid.LogGrid([c1, c2], [Int(TauGridSize / 2) + 1, TauGridSize])
    tau.grid[1], tau.grid[end] = (1.0e-8, Beta - 1.0e-8)
    return tau
end

# MomGrid construction
function kGrid(Beta, Kf, KGridSize, Type)
    @assert Type == "Fermi" || Type == "Bose" "Type should be either Fermi or Bose"
    if Type == "Fermi"
        # Fermionic Grid
        @assert MaxK > Kf "MaxK must larger than Kf!"
        kFi = floor(Int, KGridSize * log(Kf) / log(MaxK - Kf))
        lambda = sqrt(Ef * Beta) / kFi
        c1 = Grid.Coeff([0.0, Kf], [1.0, kFi + 1.0], lambda, false)
        c2 = Grid.Coeff([Kf, MaxK], [kFi, KGridSize], lambda, true)
        K = Grid.LogGrid([c1, c2], [kFi, KGridSize])
        K.grid[1] = 1.0e-6
    else
        # Bosonic Grid
        @assert MaxK > 2.0 * Kf "MaxK must larger than 2Kf!"
        kFi, twokFi = (floor(Int, KGridSize / 3), floor(Int, KGridSize / 3 * 2))
        lambda = sqrt(Ef * Beta) / kFi
        c1 = Grid.Coeff([0.0, Kf], [1.0, kFi + 1.0], lambda, true)
        c2 = Grid.Coeff([Kf, 2.0 * Kf], [kFi, twokFi + 1.0], lambda, false)
        c3 = Grid.Coeff([2.0 * Kf, MaxK], [twokFi, KGridSize], lambda, true)
        K = Grid.LogGrid([c1, c2, c3], [kFi, twokFi, KGridSize])
        K.grid[1] = 1.0e-6
    end
end

export angleGrid, tauGrid, kGrid
end
