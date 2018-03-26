Readme for Gene Experiments

In these numerical experiment we create a graph from gene expression data for the Saccharomyces cerevisiae organism and cluster genes using LambdaCC. Microarray data can be found at:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42215

We validate the clustering against Gene Ontology database (http://www.geneontology.org/page/download-annotations) and known protein-protein interections given by the String Dataset (http://version10.string-db.org/cgi/about.pl?UserId=wo13qjfRZwtl&sessionId=R5zuP5w8F9ce)

In this readme we give details for how to obtain the String annotation scores. For GO annotations we manually input clusters into the GO-slim mapper:

https://www.yeastgenome.org/cgi-bin/GO/goSlimMapper.pl

which is up to date and more complete than just the GO term data files which we pull information from.

############################################
#	Sources for Original Data	     
############################################

For all the datasets used, here I provide the URL where you can get the original information, the file the raw data was stored in, and the Matlab file I imported the data into. Usually I've just used the Matlab Import tool and extracted the relevant part of the txt file.

1. Microarray data: 
	
	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42215

	File: GSE42215_series_matrix.txt


2. Row label information can be obtained from the platform file 

	ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/GPL11232/soft/. 

	File: Platform_table.txt

	Matlab: row label data from Platform_table.txt imported into GP11232_platform_table.mat

	The main information we want is here (condensed to the most important):

	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL11232&id=12201&db=GeoDb_blob52

3. String Data: select Saccharomyces Cervisiae from dropdown menu and download interaction data file at 

	http://version10.string-db.org/cgi/download.pl?UserId=wo13qjfRZwtl&sessionId=R5zuP5w8F9ce

	File: 4932.protein.links.v10.txt


For GO annotations we simply use the GO slim mapper tool online, but many of the GO terms can be automatically downloaded:

	4. GO annotations 1: extracted from file downloaded at 

		http://www.geneontology.org/page/download-annotations
	
		File: gene_association.sgd.txt, or gene_associations_data.txt

		Matlab: go_list.mat

	5. Other source for GO associations (2): 

		http://www.yeastgenome.org/

		File: go_slim_mapping.txt

		Matlab: YeastGenome_golist.mat

	6. GO slim, parent/child relationships for GO attributes
		(this includes information about which attributes are subcategories of other attributes)

		http://geneontology.org/ontology/subsets/goslim_yeast.obo

		File: goslim_yeast.obo 


#######################################################
# Downloading, cleaning, and storing data in Matlab   
#######################################################

More information and comments are also given in Matlab script file for each step

1. Get_Microarray_Data.m 

Reading in and Cleaning the Microarray Data

a. How to read in GEO information in Matlab: 

	https://www.mathworks.com/help/bioinfo/examples/working-with-geo-series-data.html

	This gives a large matrix of microarray data—all genes are repeated more than once and some rows correspond to the ’control’ 	group, and are not genes. 

b. Load information from platform file to get row label information

c. Many rows do not correspond to genes, and all genes are repeated at least twice. Remove all rows from the ‘control’ group and for each gene in the dataset grab the first row corresponding to the gene and ignore the rest. 

d. Since each gene is repeated twice (and only 12 are repeated 4 times), we form a second matrix with the second row corresponding to each gene, just so we can compare between the two if we desire.

e. The result is two 6169 x 200 matrice G and G2 where rows are unique genes and columns are expression values for genes. Different row labels are stored in sgdid, gene_orf, and gene_symbol. Data store in GetMicroData.mat

2. Get_GO_attributes.m and Get_GO_attributes_source2.m

Reads and stores information about which Saccharomyces cerevisiae genes correspond to which GO attributes. We remove all gene-attribute pairs with the qualifier ‘NOT’, and store the rest in a (number of genes) x (number of GO_attributes) list.

There are two places we get gene-GO attribute information:

http://www.geneontology.org/page/download-annotations, File: gene_association.sgd
and

http://www.yeastgenome.org/, where under Function->Gene Ontology -> Go Slim Mapping File we are able to download the file go_slim_mapping.txt.

The second one only has 167 associated gene functions, which is very small.

3. Merge_GO_lists.m

Get_GO_attributes and Get_GO_attributes_source2 result in two gene-to-attributes dataset that only partially overlap. The lists are combined into a larger and more complete gene-to-attributes indicator matrix in this file.

4. Get_GO_comprehensive.m

The previous GO datasets don’t include information about which GO attributes are subcategories of other attributes. This .m file takes the output from Merge_GO_lists.m, GOdataMerged.mat, and fills in the gene-to-attribute matrix to include information from goslim_yeast.obo.txt regarding which GO attributes are subcategories of other GO attributes.

5. Get_String_Data

Read in gene interaction scores from String database and form a 6418 x 6418 square symmetric matrix giving scores from 0 to 1000 between different genes. Gene labels included.

#######################################################
# Running Experiments
#######################################################

1. RunGeneExperiments.m is a function that takes in just a few parameters (which half of the gene data we're using, lambda, threshold for correlation coefficients, and number of times to run lamCC algorithm)

2. Clustering_experiments.m is a script to set up parameters, call RunGeneExperiments, and save the output under a specific naming convention. For example, if we use the first half of the gene dataset (i.e. for each gene in the full microarray data, grab the first row corresponding to it), and use lambda = .99 and a threshold of .9, we cluster the data using these parameters and save it in file:

G1_Thr90_lam99.mat

3. Concise_Cluster_Analysis.m takes the output from running a clustering experiment and compares the clusters formed against GO attribute data and String data.
