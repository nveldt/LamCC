include("/homes/lveldt/GitHubRepos/corrclus/julia/grow_cliques.jl")
include("/homes/lveldt/GitHubRepos/corrclus/julia/efficient_prune_cliques.jl")
#mat = matread("../data/Rice31.mat")
using MAT
#function MainTwoAlgs(file,path,trials)
# Loading a network

trials = 1;
path = ""
file = "A_bter_1000_exp"
mat = matread(path*file*".mat")

A = mat["A"]

@printf "Starting \n"

# Two-approximation

        tic()
        Twocliques = EfficientPruneCliques.prune_cliques(A)
        #catch "Issue with 2-approx for network "*file
        #end
        TwoApptime = toc();
	n = size(A,1)
        c = ones(n,1)

        for i = 1:size(Twocliques,1)
        for j in Twocliques[i]
        c[j] = i
        end
        end

        matwrite("c_bter_1000_2app.mat",Dict( "c" => c))
