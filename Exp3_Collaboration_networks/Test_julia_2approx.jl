include("../Algs/julia/efficient_prune_cliques.jl")

using MAT

# Set number of cliques to form before grabbing the best
trials = 1;
# path = ""
# file = "ActorsMovies50"
# #or:
# #file = "AuthorPaper50"
# mat = matread(path*file*".mat")

# Load a test matrix
file = "karate"
mat = matread("../Networks/KarateA.mat")
A = mat["A"]

# Run the greedy Two-approximation (requires installation of PMC library and MaximalCliques library)
#
#
# 1. Clone PMC and MaximalCliques directories into the ../Algs folder
# 2. Make sure to compile PMC library!

        tic()
        Twocliques = EfficientPruneCliques.prune_cliques(A)
        #catch "Issue with 2-approx for network "*file
        #end
        TwoApptime = toc();
        n = size(A,1)
        b = ones(n,1)

        for i = 1:size(Twocliques,1)
        for j in Twocliques[i]
        b[j] = i
        end
        end

        open(string(file*"_"*"2approx"*".txt"),"w") do f
        write(f,file*"\n")
        write(f,"$n \n")
        write(f,"Greedy 2-approximation\n")
        write(f,"No special details\n")
        write(f,"$TwoApptime \n")
        for i = 1:n
        a = trunc(Int,b[i])
        write(f,"$a \n")
        end
        end
