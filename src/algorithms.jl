
"""
    WignerSeitz.projection(u::AbstractVector, v::AbstractVector) -> Number

Calculates the projection of `u` onto `v`: `dot(u,v) / dot(v,v)`.
"""
projection(u::AbstractVector, v::AbstractVector) = dot(u,v) / dot(v,v)

"""
    WignerSeitz.projection(M::AbstractMatrix, a, b) -> Number

Calculates the projection of `M[:,a]` onto `M[:,b]`.
"""
projection(M::AbstractMatrix, a, b) = projection(M[:,a], M[:,b])

"""
    WignerSeitz.gram(M::AbstractMatrix) -> typeof(M)

Returns the Gram matrix of M.
"""
gram(M::AbstractMatrix) = M' * M

"""
    WignerSeitz.superbase_vector(M::AbstractMatrix) -> AbstractVector

Returns the superbase vector of basis `M`, the negative sum of all basis vector (columns) combined: 
`-sum(M[:,n] for n in axes(M,2))`.
"""
superbase_vector(M::AbstractMatrix) = -sum(M[:,n] for n in axes(M,2))

"""
   WignerSeitz.superbase(M::AbstractMatrix) -> AbstractMatrix

Constructs the superbase of `M`, which is equal to `M` with an additional column containing the
negative sum of all basis vector (columns) combined: `-sum(M[:,n] for n in axes(M,2))`.
"""
superbase(M::AbstractMatrix) = hcat(M, superbase_vector(M))

"""
    WignerSeitz.selling(M::AbstractMatrix) -> typeof(M)

Returns the Selling matrix of M, which is equal to the negative Gram matrix, but includes an extra
column for the superbase vector.
"""
selling(M::AbstractMatrix) = -superbase(M)' * superbase(M)

"""
    WignerSeitz.isobtuse(M::AbstractMatrix) -> Bool

Returns `true` if a lattice basis is obtuse, meaning that all of the Selling parameters of the 
lattice basis are positive.
"""
function isobtuse(M::AbstractMatrix)
    S = superbase(M)
    for a in axes(S,2)
        for b in a:size(S,2)
            dot(S[:,a], S[:,b]) > 0 && return false
        end
    end
    return true
end
