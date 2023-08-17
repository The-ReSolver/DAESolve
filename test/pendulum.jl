# This file contains the definition of the DAE functions required to solve a
# simple pendulum motion.

# evolution and constraint functions for the DAE of a pendulum
function init_F(g, c)
    function F(q)
        x = q[1]
        y = q[2]
        u = q[3]
        v = q[4]
        λ = q[5]

        out = zeros(4)

        out[1] = u
        out[2] = v
        out[3] = -c*u - λ*x
        out[4] = -c*v - λ*y - g

        return out
    end
end
G(q) = [q[1]^2 + q[2]^2 - 1]

# ode evolution for a pendulum
function pendulum_evolution!(dθ, θ, p, t)
    # unpack parameters
    g = p[1]
    c = p[2]

    # compute rate of change of state
    dθ[1] = θ[2]
    dθ[2] = -g*sin(θ[1]) - c*θ[2]
end

# solve the pendulum using standard ODE methods
function solve_pendulum_ODE(θ₀, T; g::Real=9.81, c::Real=0.0, alg=Tsit5())
    # initialise ODE problem
    prob = ODEProblem(pendulum_evolution!, θ₀, (0, T), [g, c])

    # solve the problem
    sol = solve(prob, alg, abstol=1e-8, reltol=1e-8)
    
    return sol
end;
