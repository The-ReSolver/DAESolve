# This file defines a simple interface for a mass matrix to ensure efficient
# use of memory.

struct MassMatrix{N} <: AbstractMatrix{Int64}
    D::Vector{Int64}

    MassMatrix(D::AbstractVector) = new{length(D)}(D)
end

# default constructor for simple case where first N elements are being evolved
MassMatrix(N::Int64) = MassMatrix([n for n in 1:N])

Base.IndexStyle(::Type{<:MassMatrix}) = IndexLinear()
Base.size(::MassMatrix{N}) where {N} = (N, N)
Base.length(::MassMatrix{N}) where {N} = N
Base.getindex(M::MassMatrix, i::Int) = M.D[i]

# efficiently pick the indices of the state that are evolving
function LinearAlgebra.mul!(p::V, M::MassMatrix, q::V) where {V}
    for i in eachindex(M.D)
        p[i] = q[M.D[i]]
    end

    return p
end
