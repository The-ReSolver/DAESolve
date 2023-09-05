@testset "DAE Solver        " begin
    # setup the problem
    c = 0.1
    g = 9.81
    Δτ = 1e-4
    T = 0.2
    q1_init = Float64.([1, 0, 0, 0])
    q2_init = Float64.([0])
    θ₀ = [π/2, 0]
    F = init_F(g, c)

    # solve the DAE
    q1, q2 = solvedae(q1_init, q2_init, F, G, Δτ, maxiter=floor(Int, T/Δτ), rootopts=Options(verbose=false, use_hookstep=false))

    # solve the pendulum using ODE solvers
    sol = solve_pendulum_ODE(θ₀, T, c=c)

    # convert ODE output for comparison
    q_ode = zeros(4)
    q_ode[1] = sin(sol.u[end][1])
    q_ode[2] = -cos(sol.u[end][1])
    q_ode[3] = -q_ode[2]*sol.u[end][2]
    q_ode[4] = q_ode[1]*sol.u[end][2]

    @test q1 ≈ q_ode atol=1e-3
end
