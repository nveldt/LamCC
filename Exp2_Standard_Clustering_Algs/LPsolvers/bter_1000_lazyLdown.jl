
using MAT
include("JuliaLazyConstraints/TriFixingAlgorithms/TriFixSmallMemory.jl")
include("JuliaLazyConstraints/LamCC_LP.jl")

# A is a generated a bter graph
mat = matread("A_bter_1000_exp.mat")
A = mat["A"]

Lams = [logspace(-4,-1,10)' (.15:.1:.95)']
n = size(A,1)

open("bter1000bounds_lazyLdown.txt","w") do f
        write(f,"Clusterings for 1000 node bter graph for PNAS\n")
end

for i = 19:-1:1
    lam = Lams[i]
    D,c, lccbound = TimeFastLPlamCC(A,lam)

    println("Done with Lambda = ", lam)
    lamshort = round(lam,5)
    matwrite("bterData/D_btr_lazyLdown"*"_$lamshort.mat",Dict( "D" => D))
    matwrite("bterData/c_btr_lazyLdown"*"_$lamshort.mat",Dict( "c" => c, "lccbound" => lccbound))
    open("bter1000bounds_lazyLdown.txt","a") do f
        write(f, "$lccbound ")
        write(f, "$lamshort ")
        write(f,"\n")
    end
end
