module DAESolve

using LinearAlgebra

using NKRoots

export solvedae!

# update the state given the evolution and constraints
function _update!(q_new, q_old, F, G, Δτ, M, rootopts)
    # define an objective function from F and G
    function objective!(out, q_new)
        @views out[1:size(M, 1)] .= M*(q_new .- q_old) .+ Δτ.*F(q_new)
        @views out[(size(M, 1) + 1):end] .= G(q_new)

        return out
    end

    # perform the root finding with nkroot!
    # FIXME: redo nkroots so that it can use a pre-defined output for the objective instead of initialising a new one
    nkroot!(objective!, q_new, rootopts)

    return q_new
end

# evolve the system until a condition is met
function solvedae!(q, q_init, F, G, Δτ, M::AbstractMatrix; stopcrit=nothing, rootopts=Options(), verbose::Bool=false)
    # initialise variable for old state
    # if mass matrix isn't provided construct our own
    # loop over time steps
        # check if stopping criteria has been met
        # update the state at each time step
        # update old state with new state and go to beginning of loop
        # do some printing
    # return the solution
end
solve(q_init, F, G, Δτ; stopcrit=nothing, rootopts=Options(), verbose::Bool=false) = solve!(copy(q_init), q_init, F, G, Δτ, stopcrit=stopcrit, rootopts=rootopts, verbose=verbose)

end
