@testset "DAE Solver        " begin
    # setup the problem
    Δτ = 1e-1
    T = 20
    q_init = Float64.([1, 0, 0, 0, 0])
    F = init_F(9.81, 0.5)

    # solve the DAE
    q = solvedae(q_init, F, G, Δτ, T)

    @show q
end