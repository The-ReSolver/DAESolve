@testset "Mass Matrix       " begin
    # initialise state vectors
    q = rand(rand(5:20))
    p = zeros(rand(3:(length(q) - 1)))

    # initialise mass matrix
    M1 = MassMatrix(length(p))
    M2 = MassMatrix([rand(1:length(q)) for _ in 1:length(p)])

    mul!(p, M1, q)

    @test p == q[1:length(p)]

    mul!(p, M2, q)
    @show M2.D round.(q; digits=5) round.(p; digits=5)

    i = 1
    for j in M2.D
        @test p[i] == q[j]
        i += 1
    end
end
