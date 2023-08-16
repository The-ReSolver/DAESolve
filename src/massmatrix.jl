# This file defines a simple interface for a mass matrix to ensure efficient
# use of memory.

struct MassMatrix <: AbstractMatrix{Int64}
    D::Vector{Int64}
end

# default constructor for simple case where first N elements are being evolved
MassMatrix(N::Int64) = MassMatrix([n for n in 1:N])

# efficiently pick the indices of the state that are evolving
function LinearAlgebra.mul!(p::V, M::MassMatrix, q::V) where {V}
    for i in eachindex(M.D)
        p[i] = q[M.D[i]]
    end

    return p
end
