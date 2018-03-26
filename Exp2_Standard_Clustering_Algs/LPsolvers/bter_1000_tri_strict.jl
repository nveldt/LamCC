
using MAT
include("JuliaLazyConstraints/TriFixingAlgorithms/TriFixSmallMemory.jl")
include("JuliaLazyConstraints/LamCC_LP.jl")

# A = generate_LFR(n,15,50,.1,20,50)
mat = matread("A_bter_1000_exp.mat")
A = mat["A"]

Lams = logspace(-4,-1,10)
n = size(A,1)

open("bter_1000_trifix_bounds_strict.txt","w") do f
        write(f,"Clusterings for bter_1000 for PNAS, Tri Fix Algorithm\n")
end

Dgraph = zeros(Int8,n,n)
for i = 1:n-1
    for j = i+1:n
        if A[i,j] < .1
            Dgraph[j,i] = 1
        end
    end
end

TriTol = .1
Etol = .2

numedges = countnz(A)/2

for i = 10:-1:1
    lam = Lams[i]


    Gams = [5/lam 5/lam 10/lam 10/lam 10/lam 12/lam 12/lam 12/lam 12/lam 15/lam]

    gam = round(Gams[i],3)


    W = lam*ones(n,n)

    for i = 1:n-1
        for j = i+1:n
            if A[i,j] > .1
                W[j,i] = 1-lam
            end
        end
    end
    lamshort = round(lam,5)

    E = zeros(n,n);
    F = -gam*W;
    P = zeros(n,n)
    Q = zeros(n,n)
    M = TriangleFix(A,Dgraph,E,F,W,P,Q,Float64(gam),Float64(Etol),Float64(TriTol),lam,"TriFix_strict_$lamshort\_$gam.txt",gam)
    D = M+M'
    lccbound = sum((A[i,j]-lam)*D[i,j] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)

    println("Done with Lambda = ", lam)

    matwrite("Dbter_1000_tristrict"*"_$lamshort.mat",Dict( "D" => D))
    matwrite("bter_1000_trisctrict"*"_$lamshort.mat",Dict( "lccbound" => lccbound))
    open("bter_1000_trifix_bounds_strict.txt","a") do f
        write(f, "$lccbound ")
        write(f, "$lamshort ")
        write(f,"\n")
    end
end
