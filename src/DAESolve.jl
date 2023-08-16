module DAESolve

using LinearAlgebra

using NKRoots

export MassMatrix, solvedae!, solvedae, Options

include("massmatrix.jl")

# update the state given the evolution and constraints
# TODO: IN-PLACE!!!
function _update!(q_new, q_old, F, G, Δτ::Real, M::AbstractMatrix, rootopts::Options)
    # define an objective function from F and G
    function objective!(out, q_new)
        out[1:length(M)] .= mul!(zeros(eltype(q_new), length(M)), M, q_new .- q_old) .- Δτ.*F(q_new)
        out[(length(M) + 1):end] .= G(q_new)

        return out
    end

    # perform the root finding with nkroot!
    nkroot!(objective!, q_new, rootopts)

    return q_new
end

# evolve the system until a condition is met
function solvedae!(q, q_init, F, G, Δτ::Real, T::Real, M::Union{AbstractMatrix, Nothing}=nothing; stopcrit=(q,i)->false, rootopts::Options=Options(verbose=false), verbose::Bool=false)
    # initialise variable for old state
    q_new = copy(q_init)
    q_old = copy(q_init)

    # if mass matrix isn't provided construct our own
    if isnothing(M)
        M = MassMatrix(length(F(q_init)))
    end

    # print header if verbose
    verbose && display_header()

    # loop over time steps
    i = 0
    while i*Δτ <= T
        # update the state at each time step
        _update!(q_new, q_old, F, G, Δτ, M, rootopts)

        # check if stopping criteria has been met
        stopcrit(q_new, i) && break

        # update old state with new state and go to beginning of loop
        q_old .= q_new

        # do some printing
        verbose && display_state(q_new, F(q_new), G(q_new), i, Δτ)

        i += 1
    end
    q .= q_new

    return q
end
solvedae(q_init, F, G, Δτ, T, M=nothing; stopcrit=(q,i)->false, rootopts=Options(verbose=false, gmres_verbose=false), verbose::Bool=false) = solvedae!(similar(q_init), q_init, F, G, Δτ, T, M, stopcrit=stopcrit, rootopts=rootopts, verbose=verbose)

end
