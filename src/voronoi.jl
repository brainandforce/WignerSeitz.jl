
"""
    superbase(M::AbstractMatrix{<:Real})

Generates a superbase of `M`: a minimal basis that makes the off-diagonal elements of the Gram
matrix negative (indicating that the angles between basis vectors are maximized).
"""
function superbase(M::AbstractMatrix{<:Real})
    
end

"""
    wigner_seitz(M::AbstractMatrix{<:Real})

Generates the Wigner-Seitz cell produced by the basis vectors of a lattice.
"""
function wigner_seitz(M::AbstractMatrix{<:Real})

end
