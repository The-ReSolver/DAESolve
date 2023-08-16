# This file contains the definition of the DAE functions required to solve a
# simple pendulum motion.

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
