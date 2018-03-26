## Include software
const libpmc = joinpath(dirname(@Base.__FILE__),"libpmc")
include(joinpath(dirname(@Base.__FILE__),"..","MaximalCliques","MaximalCliques.jl"))
include(joinpath(dirname(@Base.__FILE__),"..","pmc","pmc.jl"))
##
module EfficientPruneCliques

import MaximalCliques
import PMC
using DataStructures

""" Given a set of cliques C, we repeatedly
prune the biggest clique C. And then remove those elements
from the sets.

clique_order( C::Vector{Set{Int}} ) -> Cp :: Vector{Set{Int}}

clique_order( C; minsize::Int = 1 )

Note that if minsize > 1, this will stop processing once
  the set size < minsize. This can be useful if you want
  to find a lot of cliques of size k and remove them.
  
This returns a Vector of sets. 
 """
function clique_order(C::Vector{Set{Int}}; minsize = 1)
  S = DefaultDict{Int64,Int64}(0)
  pq = PriorityQueue(Int64,Int64, Base.Order.Reverse)
  sets = DefaultDict{Int, Vector{Int}}(Vector{Int})

  Cs = Vector{Set{Int}}()

  for i in eachindex(C)
    if length(C[i]) >= minsize
      pq[i] = length(C[i])
      for j in C[i]
        push!(sets[j], i)
      end
    end
  end

  if length(pq) == 0
    return Cs
  end

  cluster = 0
  while peek(pq)[2] >= minsize
    biggest, nelem = peek(pq)
    cluster += 1 # smallest cluster is 1
    #@printf("Found clique %5i of size %4i (setid=%6i)\n", cluster, nelem, biggest)
    # remove all nodes involving biggest
    verts = Set{Int}()
    for j in C[biggest] # j is a vertex
      if S[j] == 0
        push!(verts, j)
        S[j] = cluster
        nelem -= 1
        # remove j from the sets
        for i in sets[j] # i is a set that contains vertex j
          pq[i] -= 1 # we are removing j
        end
      end
    end
    @assert nelem == 0 # we should have removed everything!
    push!(Cs, verts)
  end
  return Cs
end

function prune_cliques(A::SparseMatrixCSC;
  maximal_cliques_trial::Dates.Period = Dates.Hour(1),
  max_maximal_cliques = 100_000_000)

  n = size(A,1)
  ids = collect(1:n)
  B = copy(A)
  Cs = Vector{Set{Int}}(0)

  t0 = now() # get the current time
  t1 = t0

  maxclique = n

  while nnz(B) > 0 && maxclique > 2

    if now() - t1 >= maximal_cliques_trial
      print("running maximal cliques after ")
      print(now() - t1)
      print(" >= maximal_cliques_trial = ")
      println(maximal_cliques_trial)
      
      done, B_MCs = MaximalCliques.maximal_cliques(B; maxnum = max_maximal_cliques)
      if done
        # convert them
        println("-> maximal cliques succeeded with $(length(B_MCs)) maximal cliques")
        
        print("pruning maximal cliques")
        B_Cs = clique_order(B_MCs) # recursively prune
        for C in B_Cs
          push!(Cs, Set(map(x -> ids[x], C)))
        end
        println(" ... done")
        break
      else
        # we can still process anything equal to the current maxclique size
        println("-  maximal cliques returned $(length(B_MCs)) maximal cliques")
        B_Cs = clique_order(B_MCs; minsize = maxclique) # recursively prune anything
          # of at least the current maxclique size because you can't increase
          # clique size
        verts = [] # these are all the vertices we are going to remove
        for C in B_Cs
          push!(Cs, Set(map(x -> ids[x], C)))
          push!(verts, C...) # append all the vertices
        end
         # remove them
        left = setdiff(1:size(B,1),verts)
        ids = ids[left]
        B = B[left,left]
      end
      t1 = now()
    end

    C = PMC.pmc(B)
    maxclique = min(maxclique, length(C))

    push!(Cs, Set(map(x -> ids[x], C)))
    # MaximalCliques.maximal_cliques

    left = setdiff(1:size(B,1),C)
    ids = ids[left]
    B = B[left,left]
  end
  
  println("handling singletons and edges")

  # handle the edges and singletons
  for (ei,ej) in zip(findnz(A)...)
    push!(Cs, Set([ei,ej]))
  end
  for i=1:size(A,1)
    push!(Cs, Set([i]))
  end
  
  println("done with prune cliques")
  print("elapsed time: ", now() - t0, "\n")

  return clique_order(Cs) # handle all the singletons and edges

end

end
