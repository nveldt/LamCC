# LamCC

Code and experimental results for paper:

A Correlation Clustering Framework for Community Detection
Nate Veldt, David Gleich, Anthony Wirth

In proceedings of the 27th WWW Conference, April 2018.

[https://dl.acm.org/citation.cfm?id=3178876.3186110](https://dl.acm.org/citation.cfm?id=3178876.3186110)

Full version available at 

[https://arxiv.org/abs/1712.05825](https://arxiv.org/abs/1712.05825)

(ArXiv title: Unifying Sparsest Cut, Cluster Deletion, and Modularity Clustering Objectives with Correlation Clustering)

## Running our algorithms

Example.m is a simple file containing some examples for how to run our code on small networks. 

Additionally, we include a folder with code for each of the experiments in the paper. These experiments rely on the installation of a large number of outside packages not included in this repository. These codes may not work for you without a significant amount of setup work.

## Gurobi Software 

Note that many of our experiments use Gurobi optimization software. A free academic liscence can be obtained at Gurobi.com

## Other software

We also make use of the following, included already in our Algs folder:

* For computing adjusted rand index: [https://www.mathworks.com/matlabcentral/fileexchange/59123-density-ratio-based-clustering?focused=6842563&tab=function](http://https://www.mathworks.com/matlabcentral/fileexchange/59123-density-ratio-based-clustering?focused=6842563&tab=function)
* Setting Figure sizes: set_figure_size.m: From David F. Gleich's mcode software package. [(https://github.com/dgleich/mcode)](http://(https://github.com/dgleich/mcode))
* For cropping figures: Process_atendHeader(filein,fileout) [https://www.mathworks.com/matlabcentral/newsreader/view_thread/337705](https://www.mathworks.com/matlabcentral/newsreader/view_thread/337705)


## Experiment 5.1: Small Networks

We use a number of small datasets which can be downloaded from the SuiteSparse Matrix Collection 

[https://sparse.tamu.edu/](https://sparse.tamu.edu/)

These experiments requires code for the energy-minimization algorithm ICM of Bagon and Gallum (see [arxiv.org/abs/1112.2903](http://arxiv.org/abs/1112.2903)), available at

[github.com/shaibagon/large_scale_cc](http://github.com/shaibagon/large_scale_cc)

## Experiment 5.2: Standard Clustering Algorithms

This experiment makes use of several outside code packages. We give a brief description here. In order to obtain LP relaxations for LambdaCC, we must make use of a number of sophisticated techniques. Even then, it takes roughly a week of computation to obtain convergence.

1. **Graclus**. Graclus is a multilevel graph partitioning algorithm developed by Dhillon et al. for optimizing objectives such as normalized cut and ratio cut without computing eigenvectors. Code can be obtained at [http://www.cs.utexas.edu/users/dml/Software/graclus.html](http://www.cs.utexas.edu/users/dml/Software/graclus.html). This algorithm requires the user to specify the exact number of clusters to form. In the main text we show results for partitioning the graph into just two clusters, since this gives very good LambdaCC results for small λ. We note that nearly identical results can be obtained using the graph partitioning tool Metis [http://glaros.dtc.umn.edu/gkhome/metis/metis/overview](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).2. **Louvain**. For the local moving "Louvain" algorithm of Blondel et al., we use the implementation available at [https://github.com/ayanonagon/mkse312_project/tree/master/Community_BGLL_Matlab](https://github.com/ayanonagon/mkse312_project/tree/master/Community_BGLL_Matlab).
3. **InfoMap**. The InfoMap algorithm is similar to Louvain, but instead optimizes the map equation (see [http://www.mapequation.org/](http://www.mapequation.org/)). Source code is available at http://www.mapequation.org/code.html.4. **Recursive Maximum Quasi-Clique (RMQC)**. We employ the Quick algorithm of Liu and Wong for finding maximum quasi-cliques; their code is available at [https://www.comp.nus.edu.sg/~wongls/projects/pattern-spaces/quick-v1/](https://www.comp.nus.edu.sg/~wongls/projects/pattern-spaces/quick-v1/). We run this procedure for a specified density ρ, and extract the largest clique in the network. We then remove the quasi-clique and recurse on the remaining nodes.5. **Recursive Maximum Clique (RMC)**. For recursively extracting maximum cliques we use the Parallel Maximum Clique (PMC) library ([http://maximumclique.com/](http://maximumclique.com/)) and in some cases also call the Maximal Cliques library as a subroutine ([https://github.com/aaronmcdaid/MaximalCliques](https://github.com/aaronmcdaid/MaximalCliques))


For Lambda-Louvain, we employ the GenLouvain algorithm:

[http://netwiki.amath.unc.edu/GenLouvain/GenLouvain](http://netwiki.amath.unc.edu/GenLouvain/GenLouvain)

## Experiment 5.3: Cliques in Large Collaboration Networks

This experiment makes use of the PMC library ([http://maximumclique.com/](http://maximumclique.com/)), and the MaximalCliques library ([https://github.com/aaronmcdaid/MaximalCliques](https://github.com/aaronmcdaid/MaximalCliques)). We call both pacakges as subroutine in the Julia programming language for the recursive maximum clique algorithm.

## Experiment 5.4: Clustering Yeast Genes

Details for where to obtain the original data are included in the main text of the paper.

A number of the original data files are not included in the repository due to their size.

## Experiment 5.5: Social Network Analysis

Running Experiment 5.6 requires access to the Facebook 100 datasets, not included in this repository.

## Experiment 5.6: Clustering an Email Network

We use the Email dataset available on the SNAP repository:

[https://snap.stanford.edu/data/email-Eu-core.html](https://snap.stanford.edu/data/email-Eu-core.html)

Before running experiments we first make the edges undirected.



