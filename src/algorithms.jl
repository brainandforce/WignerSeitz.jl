
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
