module DAESolve

using LinearAlgebra, Printf

using NKRoots

export solvedae!, solvedae, Options

# TODO: make solver status a type for easier handling
# TODO: restart functionality

# update the state given the evolution and constraints
# TODO: IN-PLACE!!!
function _update!(q1, q2, q1_old, F, G, Δτ::Real, rootopts::Options)
    # define an objective function from F and G
    function objective!(out, q)
        # extract states
        q1_obj = @view(q[1:length(q1)])
        q2_obj = @view(q[(length(q1) + 1):end])

        # compute output of objective
        out[1:length(q1)] .= q1_obj .- q1_old .- Δτ.*F(q1_obj, q2_obj)
        # ? why does negating the constraints make GMRES converge much better
        out[(length(q1) + 1):end] .= .-G(q1_obj, q2_obj)

        return out
    end

    # initialise concatenated state vector
    q = vcat(q1, q2)

    # perform the root finding with nkroot!
    _, nit, gmres_maxit, gmres_minres = nkroot!(objective!, q, rootopts)

    # assign result to original state
    q1 .= q[1:length(q1)]
    q2 .= q[(length(q1) + 1):end]

    return q1, q2, nit, gmres_maxit, gmres_minres
end

# evolve the system until a condition is met
function solvedae!(q1, q2, q1_init, q2_init, F, G, Δτ::Real; maxiter::Int=1000, stopcrit=(q1, q2, i)->false, rootopts::Options=Options(verbose=false), verbose::Bool=false)
    # initialise variable for old state
    q1_old = copy(q1_init)
    q2_old = copy(q2_init)

    # set the initial state
    q1 .= q1_init
    q2 .= q2_init

    # print header if verbose
    verbose && display_header()
    verbose && display_state(0, q1, q2, F(q1, q2), G(q1, q2), NaN, NaN, 1.0)

    # loop over time steps
    i = 0
    while !(stopcrit(q1, q2, i)) && i < maxiter
        # update the state at each time step
        _, _, nit, gmres_maxit, gmres_minres = _update!(q1, q2, q1_old, F, G, Δτ, rootopts)

        # update old state with new state and go to beginning of loop
        q1_old .= q1
        q2_old .= q2

        # update iterator
        i += 1

        # do some printing
        verbose && display_state(i, q1, q2, F(q1, q2), G(q1, q2), nit, gmres_maxit, gmres_minres)
    end

    return q1, q2
end
solvedae(q1_init, q2_init, F, G, Δτ; maxiter::Int=1000, stopcrit=(q1, q2, i)->false, rootopts=Options(verbose=false, gmres_verbose=false), verbose::Bool=false) = solvedae!(similar(q1_init), similar(q2_init), q1_init, q2_init, F, G, Δτ, maxiter=maxiter, stopcrit=stopcrit, rootopts=rootopts, verbose=verbose)

function display_header()
    println("----------------------------------------------------------------------------------------------------------------")
    println("|   iter  |    ||F||    |    ||G||    |    ||r||    |   Newton Iter   |    Max GMRES Iter   |   Min GMRES Res  |")
    println("----------------------------------------------------------------------------------------------------------------")
    flush(stdout)
    return nothing
end

function display_state(i, q1, q2, F, G, nit, git, grnorm)
    str = @sprintf("|  %4d   |   %5.2e  |   %5.2e  |   %5.2e  |     %4d        |       %4d          |     %5.2e     |", i, norm(F), norm(G), norm(@view(q2[1:3])), nit, git, grnorm)
    println(str)
    flush(stdout)
    return nothing
end

end
