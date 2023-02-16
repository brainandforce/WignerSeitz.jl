module WignerSeitz

using LinearAlgebra
using Requires
import NormalForms: eye, detb!, detb, isunimodular

include("algorithms.jl")
include("reduction.jl")
export LatticeReduction, HermiteReducedLattice
export hermite_reduce!, hermite_reduce

#=
include("voronoi.jl")
export superbase, wigner_seitz
=#

function __init__()
    # Contains methods specific to SMatrix, which can't be mutated in-place
    @require StaticArrays="90137ffa-7385-5640-81b9-e52037218182" begin
        include("staticarrays.jl")
    end
end

end
