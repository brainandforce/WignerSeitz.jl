using .StaticArrays

function hermite_reduce(M::SMatrix{D1,D2,T}) where{D1,D2,T}
    return SMatrix{D1,D2,T}(hermite_reduce!(convert(Matrix, M)))
end
