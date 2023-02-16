
"""
    LatticeReduction{T,M<:AbstractMatrix{T}} <: LinearAlgebra.Factorization{T}

Supertype for the results of lattice reductions: algorithms which convert an `m`×`n` matrix (which
represents an `n`-dimensional lattice in an `m`-dimensional space) to a sufficiently short basis.
"""
abstract type LatticeReduction{T,M<:AbstractMatrix{T}} <: LinearAlgebra.Factorization{T}
end

# Convert the unimodular factor to an integer type
function Base.getproperty(x::LatticeReduction, s::Symbol)
    return s === :U ? Int.(getfield(x, :U)) : getfield(x, s)
end

LinearAlgebra.issuccess(H::LatticeReduction) = iszero(H.info)

Base.iterate(L::LatticeReduction) = (L.B, Val{:U}())
Base.iterate(L::LatticeReduction, ::Val{:H}) = iterate(L)
Base.iterate(L::LatticeReduction, ::Val{:U}) = (L.U, nothing)
Base.iterate(L::LatticeReduction, ::Nothing) = nothing

function Base.summary(io::IO, L::LatticeReduction)
    print(io, join(size(L.B), '×'), ' ', typeof(L), ":")
end

function Base.show(io::IO, mime::MIME"text/plain", L::LatticeReduction)
    if issuccess(L)
        summary(io, L)
        println(io, "\nReduced basis:")
        Base.show(io, mime, L.B)
        println(io, "\nUnimodular factor:")
        Base.show(io, mime, L.U)
    else
        print(io, "Failed factorization of type $(typeof(L))")
    end
end

"""
    HermiteReducedLattice{T,M<:AbstractMatrix{T}} <: LatticeReduction{T,M}

Stores a Hermite reduced basis `B` and the associated unimodular factor `U` such that for an input
matrix `M`, `L = U*M`. This is an exactly reduced lattice.
"""
struct HermiteReducedLattice{T,M} <: LatticeReduction{T,M}
    B::M
    U::M
    info::LinearAlgebra.BlasInt
    function HermiteReducedLattice(B::AbstractMatrix, U::AbstractMatrix)
        M = Base.promote_typeof(B, U)
        return new{eltype(M),M}(B, U, (isunimodular(U) ? 0 : -1))
    end
end

"""
    hermite_reduce!(M::AbstractMatrix) -> HermiteReducedLattice{eltype{M},M}

Performs an exact lattice reduction (in other words, without approximations to perform the reduction
in polynomial time) of a lattice whose basis vectors are the columns of `M`. This lattice reduction
occurs in-place, and the result is returned as a `HermiteReducedLattice`.

For a lattice reduction which creates a copy of the original matrix, or lattice reduction of an
`SMatrix`, see `hermite_reduction`.
"""
function hermite_reduce!(M::AbstractMatrix)
    X = copy(M)
    # Right unimodular matrix factor
    U = eye(M, 2)
    # Size rank of the vector we're minimizing against
    vectors_minimized = falses(size(M,2))
    # Norms of all the vectors
    norms = [norm(M[:,n]) for n in axes(M,2)]
    # Don't stop until we minimized everything
    while !all(vectors_minimized)
        norm_order = sortperm(norms)
        # Start with the shortest vector that hasn't been minimized
        a = norm_order[count(vectors_minimized) + 1]
        for b in filter(!isequal(a), norm_order)
            # Shorten all of the other vectors and recalculate the norm
            x = round(projection(M, b, a))
            M[:,b] .-= M[:,a] * x
            U[:,b] .-= U[:,a] * x
            norms[b] = norm(M[:,b])
        end
        # Report if the current vector was minimized
        vectors_minimized[a] = sortperm(norms)[count(vectors_minimized) + 1] == a
        # Mark vectors that changed ordering as unminimized
        vectors_minimized .&= (norm_order .== sortperm(norms))
    end
    return HermiteReducedLattice(M, U)
end

"""
    hermite_reduce(M::AbstractMatrix) -> HermiteReducedLattice{eltype{M},M}

Performs an exact lattice reduction (in other words, without approximations to perform the reduction
in polynomial time) of a lattice whose basis vectors are the columns of `M`. This creates a copy of
the input matrix and modifies that copy in-place.
"""
hermite_reduce(M::AbstractMatrix) = hermite_reduce!(deepcopy(M))
