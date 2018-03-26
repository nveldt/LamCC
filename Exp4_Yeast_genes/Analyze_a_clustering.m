%% Load and analyze a clustering

load Genes131opt
n = numel(c);

fprintf('Looking at a clustering of %d nodes.\n',n);

% Set FullDetails = 1 if you want to see full details for each cluster
% (names, all GO terms associated, etc.)
FullDetails = 1;

% There are 4 ways to valiadate against GO attribute information:
% GoData = 1 means check against data from gene_association.sgd
% GoData = 2 means check againt data from go_slim_mapping.txt 
% GoData = 3 means check against data from both 1 and 2
% GoData = 4 means check against both 1 and 2, and additionally include
%            information about which GO attributes are subcategories of other GO
%            attributes
GOdata = 1;    

% Set the minimum size clique we want to consider
MinSize = 3;

% NOTE: You can get more complete GO annotation information by taking the
% sets of genes in each cluster and using the GO slim mapper tool:
%
%       https://www.yeastgenome.org/cgi-bin/GO/goSlimMapper.pl
%
% This is what we do in our paper.

%% Get the set largest cliques

numNodes = zeros(max(c),1);
clusterID = (1:max(c))';
for i = 1:max(c)
    numNodes(i) = numel(find(c == i));
end

LargeCliques = find(numNodes >= MinSize);
CliqueInfo = [LargeCliques numNodes(LargeCliques)];


%% Check the GO attributes for each cluster

switch GOdata
    case 1
        % From: gene_association.sgd
        % Processed in Get_GO_attributes
        load MatFiles/GOdata1.mat
    case 2
        % From: go_slim_mapping.txt
        % Processed in Get_GO_attributes_source2
        load MatFiles/GOdata2.mat 
    case 3
        % If we merge GOdata1 and GOdata2
        load MatFiles/GOdataMerged.mat
    otherwise
        % If we take GOdataMerged, and include information about which GO functions
        % are subcategories of other attributes, see file goslim_yeast.obo.txt
        load MatFiles/AllGOattributes.mat
        GeneToGO = GeneToFullGO;
end


%% Print out concise stats about each cluster

load MatFiles/StringMatrix
h = size(A,1);

% Insert an extra space for whenever we deal with a gene in our experiment
% that has no correponding data in the string matrix
A = [A zeros(h,1); zeros(1,h+1)];
string_labels = [string_labels; 'NoDataGene'];

fprintf('\nConcise Clustering Information\n')
fprintf('We give cluster size, number of GO attributes that all share,\n')
fprintf('percent of all genes in the dataset associated with the most rare of these attributes,\n')
fprintf('And the Min, Max, and Avg String score among all pairs of genes in the cluster\n\n')
fprintf('Clus  Size\t#GO \t %% \tMin  \t Max \t Avg \n')

for i = 1:size(CliqueInfo,1)
    ID = CliqueInfo(i,1);       % The clique ID in the clustering vector c
    Members = find(c == ID);    % The indices in c in this clique

    Local_genes = unique(SGD(Members)); % SGD gene names for this clique
    geneIDs = [];
    cliqueSize = numel(Local_genes);
    
    % 'genes' gives row/gene labels for the GeneToGO matrix
    % For each gene in this clique, find its ID in the GeneToGO matrix
    for j = 1:numel(Local_genes)
        geneIDs = [geneIDs; find(strcmp(genes,Local_genes(j)))];
    end
    geneIDs = unique(geneIDs);
    GOmat = GeneToGO(geneIDs,:);     % Indicator matrix for which GO attributes each gene in the cluster corresponds to
    GOids = find(sum(GOmat,1));      % All GO attributes associated with at least one gene in this cluster
    LocalGOlist = GO_list(GOids);    % names of the associated GO attributes
    Local_GeneToGo = GOmat(:,GOids); % Local gene-to-attribute indicator matrix
    
    GoCount = sum(Local_GeneToGo,1);           % tells me how many genes in this cluster share the attribute
    AllShare = find(GoCount == cliqueSize);    % which GO attributes are shared by all genes
    
    
    TotalGoCount = sum(GeneToGO(:,GOids),1);   % how many genes in the data we have share the attribute
    Percents = TotalGoCount./numel(genes)*100; % Percentage of the entire dataset that exhibit this GO attribute
    Percents = Percents(AllShare);             % We care only about GO attributes shared by all genes in the cluster
    
    [BestPercent,wherebest] = min(Percents);
    
    % Now we want to find string data
    % We must compare ORF names in this cluster
    % against data in the String score matrix A
    % 'genes' gives row/gene labels for the GeneToGO matrix
    % For each gene in this clique, find its ID in the GeneToGO matrix
    
    
    Local_orfs = unique(ORF(Members)); % ORF gene names for this clique
    orfIDs = [];
    for j = 1:numel(Local_orfs)
        place = find(strcmp(string_labels,Local_orfs(j)));
        if numel(place) > 0
            orfIDs = [orfIDs; place];
        else
            orfIDs = [orfIDs; h+1];
        end
    end
    orfIDs = unique(orfIDs);
    if cliqueSize < 2
        stringVec = 0;
    else
    Local_StringMat = A(orfIDs,orfIDs);
    % Now vectorize the string matrix
    stringVec = nchoosek(cliqueSize,2);
    next = 0;
    for col = 2:cliqueSize
        for row = 1:col-1
            next = next+1;
            stringVec(next) = Local_StringMat(row,col);
        end
    end
    end
    % Output data
    if cliqueSize >= MinSize
    fprintf(' %d \t',i)
    fprintf('%d \t',cliqueSize)
    fprintf('%d \t',full(nnz(AllShare)))
    fprintf('%.1f \t',full(BestPercent))
    fprintf('%d \t %d \t %.2f', min(stringVec),max(stringVec),mean(stringVec))
    fprintf('\n')
    end
end


%% Give more detailed information about each cluster
if FullDetails
fprintf('\nConsider each cluster individually in more detail. \n')
fprintf('This includes information about the GO terms for the cluster. \n')
fprintf('For more complete information about cluster GO terms, input the names of the genes into the SGD GO slim mapper tool: \n')
fprintf('https://www.yeastgenome.org/cgi-bin/GO/goSlimMapper.pl \n')


% Yes or no: show attribute if all but one genes shares it
minusOnelist = 0;  

% When printing results for GO attributes not shared by all nodes,
% only show the GO attribute if it is 'rare'
percentCut = 100;

for i = 1:size(CliqueInfo,1)
    ID = CliqueInfo(i,1);       % The clique ID in the clustering vector c
    Members = find(c == ID);    % The indices in c in this clique

    Local_genes = SGD(Members); % SGD gene names for this clique, for finding GO terms
    sys_names = ORF(Members);   % We'll print out systematic names
    geneIDs = [];
    cliqueSize = numel(Members);
    
    % 'genes' gives row/gene labels for the GeneToGO matrix
    % For each gene in this clique, find its ID in the GeneToGO matrix
    for j = 1:numel(Local_genes)
        geneIDs = [geneIDs; find(strcmp(genes,Local_genes(j)))];
    end
    
    GOmat = GeneToGO(geneIDs,:);     % Indicator matrix for which GO attributes each gene in the cluster corresponds to
    GOids = find(sum(GOmat,1));      % All GO attributes associated with at least one gene in this cluster
    LocalGOlist = GO_list(GOids);    % names of the associated GO attributes
    Local_GeneToGo = GOmat(:,GOids); % Local gene-to-attribute indicator matrix
    
    GoCount = sum(Local_GeneToGo,1);           % tells me how many genes in this cluster share the attribute
    TotalGoCount = sum(GeneToGO(:,GOids),1);   % how many genes in the data we have share the attribute
    
    % Output data
    fprintf('\nProfile for Cluster %d \n',i)
    fprintf('----------------------------\n')
    fprintf('Number of genes = %d \n',numel(Members));
    fprintf('Gene names: \n')
    for g = 1: cliqueSize;
        fprintf('%s \n',char(sys_names(g)))
    end

    fprintf('GO attribute data \n')
    for t = 1:numel(LocalGOlist)
        
       totalpercent = full(TotalGoCount(t)/numel(genes)*100);
       if GoCount(t) == cliqueSize
       fprintf('%d of %d genes share attribute %s ',full(GoCount(t)),cliqueSize,char(LocalGOlist(t)))
       fprintf(' (Present in %.3f percent of all genes)\n',totalpercent)
       end
       
       if GoCount(t) == cliqueSize -1 && minusOnelist
           fprintf('%d of %d genes share attribute %s ',full(GoCount(t)),cliqueSize,char(LocalGOlist(t)))
           fprintf(' (Present in %.3f percent of all genes)\n',totalpercent)
       end
    end
end
end