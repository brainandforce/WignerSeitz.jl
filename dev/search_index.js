var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = WignerSeitz","category":"page"},{"location":"#WignerSeitz","page":"Home","title":"WignerSeitz","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for WignerSeitz.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [WignerSeitz]","category":"page"},{"location":"#WignerSeitz.HermiteReducedLattice","page":"Home","title":"WignerSeitz.HermiteReducedLattice","text":"HermiteReducedLattice{T,M<:AbstractMatrix{T}} <: LatticeReduction{T,M}\n\nStores a Hermite reduced basis B and the associated unimodular factor U such that for an input matrix M, L = U*M. This is an exactly reduced lattice.\n\n\n\n\n\n","category":"type"},{"location":"#WignerSeitz.LatticeReduction","page":"Home","title":"WignerSeitz.LatticeReduction","text":"LatticeReduction{T,M<:AbstractMatrix{T}} <: LinearAlgebra.Factorization{T}\n\nSupertype for the results of lattice reductions: algorithms which convert an m×n matrix (which represents an n-dimensional lattice in an m-dimensional space) to a sufficiently short basis.\n\n\n\n\n\n","category":"type"},{"location":"#WignerSeitz.gram-Tuple{AbstractMatrix}","page":"Home","title":"WignerSeitz.gram","text":"WignerSeitz.gram(M::AbstractMatrix) -> typeof(M)\n\nReturns the Gram matrix of M.\n\n\n\n\n\n","category":"method"},{"location":"#WignerSeitz.hermite_reduce!-Tuple{AbstractMatrix}","page":"Home","title":"WignerSeitz.hermite_reduce!","text":"hermite_reduce!(M::AbstractMatrix) -> HermiteReducedLattice{eltype{M},M}\n\nPerforms an exact lattice reduction (in other words, without approximations to perform the reduction in polynomial time) of a lattice whose basis vectors are the columns of M. This lattice reduction occurs in-place, and the result is returned as a HermiteReducedLattice.\n\nFor a lattice reduction which creates a copy of the original matrix, or lattice reduction of an SMatrix, see hermite_reduction.\n\n\n\n\n\n","category":"method"},{"location":"#WignerSeitz.hermite_reduce-Tuple{AbstractMatrix}","page":"Home","title":"WignerSeitz.hermite_reduce","text":"hermite_reduce(M::AbstractMatrix) -> HermiteReducedLattice{eltype{M},M}\n\nPerforms an exact lattice reduction (in other words, without approximations to perform the reduction in polynomial time) of a lattice whose basis vectors are the columns of M. This creates a copy of the input matrix and modifies that copy in-place.\n\n\n\n\n\n","category":"method"},{"location":"#WignerSeitz.projection-Tuple{AbstractMatrix, Any, Any}","page":"Home","title":"WignerSeitz.projection","text":"WignerSeitz.projection(M::AbstractMatrix, a, b) -> Number\n\nCalculates the projection of M[:,a] onto M[:,b].\n\n\n\n\n\n","category":"method"},{"location":"#WignerSeitz.projection-Tuple{AbstractVector, AbstractVector}","page":"Home","title":"WignerSeitz.projection","text":"WignerSeitz.projection(u::AbstractVector, v::AbstractVector) -> Number\n\nCalculates the projection of u onto v: dot(u,v) / dot(v,v).\n\n\n\n\n\n","category":"method"},{"location":"#WignerSeitz.selling-Tuple{AbstractMatrix}","page":"Home","title":"WignerSeitz.selling","text":"WignerSeitz.selling(M::AbstractMatrix) -> typeof(M)\n\nReturns the Selling matrix of M, which is equal to the negative Gram matrix, but includes an extra column for the superbase vector.\n\n\n\n\n\n","category":"method"},{"location":"#WignerSeitz.superbase-Tuple{AbstractMatrix}","page":"Home","title":"WignerSeitz.superbase","text":"WignerSeitz.superbase(M::AbstractMatrix) -> AbstractMatrix\n\nConstructs the superbase of M, which is equal to M with an additional column containing the negative sum of all basis vector (columns) combined: -sum(M[:,n] for n in axes(M,2)).\n\n\n\n\n\n","category":"method"},{"location":"#WignerSeitz.superbase_vector-Tuple{AbstractMatrix}","page":"Home","title":"WignerSeitz.superbase_vector","text":"WignerSeitz.superbase_vector(M::AbstractMatrix) -> AbstractVector\n\nReturns the superbase vector of basis M, the negative sum of all basis vector (columns) combined:  -sum(M[:,n] for n in axes(M,2)).\n\n\n\n\n\n","category":"method"}]
}
