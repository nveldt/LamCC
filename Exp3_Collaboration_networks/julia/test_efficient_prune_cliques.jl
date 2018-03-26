include("efficient_prune_cliques.jl")
using Base.Test

@testset "Prune Clique Tests" begin

    A = ones(5,5)
    B = ones(3,3)
    C = ones(2,2)
    G = blkdiag(sparse(A),sparse(B),sparse(C))
    
    @testset "Repeated pmc calls" begin
        C = EfficientPruneCliques.prune_cliques(G)
        @test length(C[1]) == 5
        @test length(C[2]) == 3
        @test length(C[3]) == 2
        @test length(setdiff(C[1],1:5)) == 0
        @test length(setdiff(C[2],6:8)) == 0
        @test length(setdiff(C[3],9:10)) == 0
    end
    
    @testset "Maximal clique only calls" begin
        C = EfficientPruneCliques.prune_cliques(G;maximal_cliques_trial = Dates.Second(0))
        @test length(C[1]) == 5
        @test length(C[2]) == 3
        @test length(C[3]) == 2
        @test length(setdiff(C[1],1:5)) == 0
        @test length(setdiff(C[2],6:8)) == 0
        @test length(setdiff(C[3],9:10)) == 0
    end
    
    @testset "Mixed call" begin
        C = EfficientPruneCliques.prune_cliques(G;maximal_cliques_trial = Dates.Millisecond(1))
        @test length(C[1]) == 5
        @test length(C[2]) == 3
        @test length(C[3]) == 2
        @test length(setdiff(C[1],1:5)) == 0
        @test length(setdiff(C[2],6:8)) == 0
        @test length(setdiff(C[3],9:10)) == 0
    end
    
    @testset "Mixed call partial solution" begin
        A = ones(5,5)
        B = ones(5,5)
        C = ones(3,3)
        D = ones(3,3)
        G = blkdiag(sparse(A),sparse(B),sparse(C),sparse(D))
        C = EfficientPruneCliques.prune_cliques(G;
            maximal_cliques_trial = Dates.Millisecond(1),
            max_maximal_cliques = 2)
        @test length(C[1]) == 5
        @test length(C[2]) == 5
        @test length(C[3]) == 3
        @test length(C[4]) == 3
    end
    
    
end 