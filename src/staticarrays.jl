using .StaticArrays

function hermite_reduce(M::SMatrix{D1,D2}) where{D1,D2}
    L = hermite_reduce!(convert(Matrix, M))
    return HermiteReducedLattice(SMatrix{D1,D2}(L.B), SMatrix{D2,D2}(getfield(L, :U)), L.info)
end
