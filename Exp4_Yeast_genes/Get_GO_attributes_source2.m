%% Now we associate each gene with its GO functions
%
% These can be found at http://www.yeastgenome.org/, under
% Function->Gene Ontology -> Go Slim Mapping File
%
% I just want the name of the gene and the corresponding GO attributes, so
% I imported those through the Matlab Import Data tool

clear

% gives a list of gene-attribute pairs

load MatFiles/YeastGenome_golist    


%% Find out how many GO-attributes there are

FullAttList = go_list(:,3);
GO_list = unique(FullAttList);
GO_list = GO_list(2:end);           % get rid of empty string entry
numGOtags = numel(GO_list)          % There are 167 GO-attributes
                            
%% e.g. GO:0008150 is 91
k = 'GO:0008150';
where = find(strcmp(GO_list,k))
fprintf('%s is attribute number %d in the list \n',char(k),where)

%% Find out how many genes there are

FullGeneList = go_list(:,1);
genes = unique(FullGeneList);
numGenes = numel(genes)         % There are 7013 genes

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

save('GOdata2','GO_list','genes','GeneToGO')