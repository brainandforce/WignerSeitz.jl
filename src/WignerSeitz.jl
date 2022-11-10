module WignerSeitz

using LinearAlgebra

"""
    reduce_basis!(M::AbstractMatrix{<:Real})

Reduces the basis vectors associated with a lattice to a form where all lattice vectors have
minimum length and all angles between lattice vectors are at least 90 degrees. This is the form of
the "superbase" described in Conway and Sloane's 1992 paper.

Essentially, this can be thought of as a "QRU" decomposition, where `Q` is a point isometry 
(rotation/reflection/etc.), `R` is an upper triangular matrix (this essentially reorients the basis
vectors of the space), and `U` is a unimodular matrix (all lattices whose basis vectors are 
transformed by a unimodular matrix are identical), which permutes the columns of `R`, In this case,
we minimize the dot product ofeach column of `R` with itself by right multiplication with a 
suitable unimodular matrix.
"""
function reduce_basis!(M::AbstractMatrix{<:Real})
    # Perform a QR decomposition to reorient the basis vectors
    (Q,R) = qr!(M)
    # Unimodular factor starts as an identity matrix
    U = Matrix{Int}(LinearAlgebra.I(size(R,1)))
    # Reduce the off-diagonal elements of R to the smallest non-positive numbers possible
    # It's easier if we do this backwards - also, we can skip the first column
    for a in reverse(minimum(axes(R)))
        # Now minimize the off-diagonal elements - again, do this backwards
        for b in reverse(axes(R,1)[1:a-1])
            # Construct mul so we minimize abs(R[a-1,a])
            mul = Int(div(R[b,a], R[b,b], RoundNearest))
            for c in axes(R,1)
                R[c,a] -= mul * R[c,b]
                U[b,c] += mul * U[a,c]
            end
            @info "At a = $a, b = $b:" *
                "\nmul = $mul" * 
                "\nR = " * repr("text/plain", R) *
                "\nU = " * repr("text/plain", U) *
                "\nNorms: " * repr(norm.(eachcol(R)))
        end
        # Invert column so the sign of the diagonal element is positive
        s = _sign(R[a,a])
        for c in axes(R,1)
            R[c,a] *= s
            U[a,c] *= s
        end
    end
    # Perform operations to make each column have a negative dot product with all others
    # This is restricted to inversion or adding/subtracting one column to/from another

    return (Q = Q, R = R, U =U)
end

reduce_basis(M::AbstractMatrix{<:Real}) = reduce_basis!(deepcopy(M))

#=
"""
    wigner_seitz(M::AbstractMatrix{<:Real})

Finds the Wigner-Seitz cell produced by the basis vectors of a 
"""
function wigner_seitz(M::AbstractMatrix{<:Real})

end

export wigner_seitz
=#

end
