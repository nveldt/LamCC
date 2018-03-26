%% Read microarray data from file
clear
data = geoseriesread('RawData/GSE42215_series_matrix.txt');

%% Get row labels 

load MatFiles/GPL11232_platform_table

% Data imported from part of the platform file using Matlab's import tool
% 
% This is exactly the information that can be found at:
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL11232&id=12201&db=GeoDb_blob52

%  Now we want to extract the most useful information:
%  Gene expression level data matrix, with no control rows
%  Most relevent row label information and column label information
%  In particular, the 8th row label tells you whether this row is gene or
%  control.

%% Remove 'control' rows

A = data.Data*eye(200);             % Turn into just a regular matrix, not bioma.data.DataMatrix
n = size(A,1);
GeneIndicator = zeros(n,1);
ReporterGroup = row_labels(2:end,8);

for i = 1:n
    if strcmp(cellstr(ReporterGroup(i)),'gene')
        GeneIndicator(i) = 1;
    end
        
end
nnz(GeneIndicator)                  % still more rows than genes
GeneRows = find(GeneIndicator);

%% Count how many times each gene shows up

sgdid = row_labels(2:end,12);
orf = row_labels(2:end,10);

sgdid = sgdid(GeneRows);
orf = orf(GeneRows);   

%% There are two rows that have an associated systematic name
% but their SGD id is not listed. I will fix that here so we can have both
% labels for future use.
% Speficially, the gene is YGR226C, and its SGD id is S000003458,
% as can be seen here: http://www.yeastgenome.org/locus/S000003458/overview

sgdid(4328)
sgdid(4189)
orf(4328)
orf(4189)

sgdid(4328) = cellstr('S000003458');
sgdid(4189) = cellstr('S000003458');

%% Confirm that all genes have two kinds of labels

for i = 1:numel(orf)
    if strcmp(sgdid(i),'') || strcmp(orf(i),'')
        i
        sgdid(i)
        orf(i)
    end
    if numel(orf(i)) == 0
        i
    end
end

%%

m = 40000;
counter = zeros(m,1);
for i = 1:(size(sgdid,1))
    str = char(sgdid(i));
    num = str2num(str(2:end));
    counter(num) = counter(num)+1;
end
nnz(counter)     % tells you that there are 6170 unique cell ids

%% 12 genes are repeated 4 times. 6158 are repeated twice.
nnz(counter == 4)
nnz(counter == 2)
uniqueGeneNum = nnz(counter > 0)

%% Get the first and second time each gene shows up

FirstInstance = zeros(uniqueGeneNum,1);
SecondInstance = zeros(uniqueGeneNum,1);

place = 0;

thegenes = unique(sgdid);

for i = 1:uniqueGeneNum  
        place = place+1;
        str = char(thegenes(i));
        Times = find(strcmp(sgdid,str));
        
        % confirms there are always at least 2 rows for each gene
        assert(numel(Times) > 1); 
        
        FirstInstance(place) = Times(1);
        SecondInstance(place) = Times(2);
end

% we end up with 6170 unique genes, repeated twice
numel(FirstInstance) 
numel(SecondInstance)

%% Get the global indices--
% (currently we only have indices within the 12364 rows corresponding to
% genes)

UniqueIDs1 = FirstInstance;
Genes1 = GeneRows(UniqueIDs1); 

UniqueIDs2 = SecondInstance;
Genes2 = GeneRows(UniqueIDs2);

%% Labels are now unique

labels = row_labels(2:end,:);         % all label info except title
labels1 = labels(Genes1,:);           % Gene dataset 1
labels2 = labels(Genes2,:);           % Gene dataset 2

size(labels1)
size(labels2)
G1 = A(Genes1,:);
G2 = A(Genes2,:);

% These will be the row label for both matrices--the first and second time
% a gene shows up
gene_orf = orf(FirstInstance);       
gene_sgdid = sgdid(FirstInstance);

%% Save the output, two 6170 x 200 matrices with the same sgdid and orf row labels
save('MatFiles/GetMicroData.mat','G1','G2','gene_sgdid','gene_orf')