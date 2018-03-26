%% Now we associate each gene with its GO functions


% These can be found at http://www.geneontology.org/page/download-annotations
%
% I just want the name of the gene and the corresponding GO attributes, so
% I imported those through the Matlab Import Data tool

clear

% go_list.mat gives a list of gene-attribute pairs
load MatFiles/go_list.mat       
    

%% Remove genes in the list explicitely 'NOT' associated with a function

notes = (go_list(:,2));
Not = find(strcmp(notes,'NOT'));
All = 1:size(go_list,1);
GoodAssociations = setdiff(All,Not);
numel(GoodAssociations)
go_list = go_list(GoodAssociations,:);

%% Find out how many GO-attributes there are

FullAttList = go_list(:,3);
GO_list = unique(FullAttList);
numGOtags = numel(GO_list)      % There are 5767 GO-attributes
                            
%% e.g. GO:0000278 is 158 in the list, GO:0008150 is 1993
k = 'GO:0008150';
where = find(strcmp(GO_list,k))
fprintf('%s is attribute number %d in the list \n',char(k),where)

%% Find out how many genes there are

FullGeneList = go_list(:,1);
genes = unique(FullGeneList);
numGenes = numel(genes)         % There are 6448 genes

%% For every gene in the associations file, I want to extract its GO attributes
% I'll put them all in a (number of genes) x (number of attributes)
% indicator matrix, and create both row labels and column labels

GeneToGO = zeros(numGenes,numGOtags);

for i = 1:numel(go_list(:,1))
   thegene = char(go_list(i,1));
   GeneID = find(strcmp(genes,thegene));
   gotag = char(go_list(i,3));
   TagID = find(strcmp(GO_list,gotag));
   GeneToGO(GeneID,TagID) = 1;
end

%% Save the info

save('GOdata1','GO_list','genes','GeneToGO')