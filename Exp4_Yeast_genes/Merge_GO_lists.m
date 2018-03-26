%% We want to combine GeneToGO data from the two different sources

clear
load MatFiles/GOdata1.mat
golist1 = GO_list;
genes1 = genes;
GtG1 = GeneToGO;
numgenes1 = numel(genes1)
numgo1 = numel(golist1)

load MatFiles/GOdata2.mat
golist2 = GO_list;
genes2 = genes;
GtG2 = GeneToGO;
numgenes2 = numel(genes2)
numgo2 = numel(golist2)

%% If we first change SGD identifiers and GO attributes to integers, it 
% will be easier to work with

GoIntList1 = zeros(numgo1,1);

for i = 1:numgo1
    gonum = char(golist1(i));
    gnum = str2num(gonum(4:end));
    GoIntList1(i) = gnum;
end

sgdInts1 = zeros(numgenes1,1);

for i = 1:numgenes1
    genesgd = char(genes1(i));
    sgdnum = str2num(genesgd(2:end));
    sgdInts1(i) = sgdnum;
end

%% Same for the second dataset

GoIntList2 = zeros(numgo2,1);

for i = 1:numgo2
    gonum = char(golist2(i));
    gnum = str2num(gonum(4:end));
    GoIntList2(i) = gnum;
end

sgdInts2 = zeros(numgenes2,1);

for i = 1:numgenes2
    genesgd = char(genes2(i));
    sgdnum = str2num(genesgd(2:end));
    sgdInts2(i) = sgdnum;
end

%% How many genes are in common vs total?

fullGeneSet = union(sgdInts2,sgdInts1);
TotalGenes = numel(fullGeneSet)
numel(intersect(sgdInts2,sgdInts1))

fullGoSet = union(GoIntList1,GoIntList2);
TotalGo = numel(union(GoIntList1,GoIntList2))
numel(intersect(GoIntList1,GoIntList2))

%% Actually let's stick with the non-integer stuff

fullGeneSet = union(genes1,genes2);
TotalGenes = numel(fullGeneSet)

fullGoSet = union(golist1,golist2);
TotalGo = numel(fullGoSet)


%% Combine the information

GeneToGO = zeros(TotalGenes,TotalGo);

% Iterate through every gene, and get its corresponding GO attriutes from
% two places, put them in a list, and find their indices in the fullGoSet
% list. Then update the GeneToGo matrix accordingly

for j = 1:TotalGenes
    gene = fullGeneSet(j);
    
    % Get the go attributes it corresponds to from the first dataset
    GeneID1 = find(strcmp(genes1,gene));
    First = find(GtG1(GeneID1,:));
    
    GeneID2 = find(strcmp(genes2,gene));
    Second = find(GtG2(GeneID2,:));
   
    GO1 = golist1(First);
    GO2 = golist2(Second);
    
    FullSet = unique([GO1;GO2]);
    
    goids = zeros(numel(FullSet),1);
    for k = 1:numel(FullSet)
        goids(k) = find(strcmp(fullGoSet,FullSet(k)));
    end
    
    GeneToGO(j,goids) = 1;
        
end

%%
GO_list = fullGoSet;
genes = fullGeneSet;
save('GOdataMerged','GO_list','genes','GeneToGO')
