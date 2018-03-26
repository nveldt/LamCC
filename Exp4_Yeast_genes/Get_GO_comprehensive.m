%% Load the current GO attributes data
clear
load MatFiles/GOdataMerged.mat

% 'GeneToGO' is an indicator matrix for which genes have which GO attribute,
% based on two different data files (this is the output from
% Merge_GO_lists)
%
% 'GO_list' gives column labels 
% 'genes' gives row labels

%% We will build a GO attribute ancestor-descendent graph
n = numel(GO_list);
gograph = zeros(n);     % create a non-symmetric graph of relationships
                        % between parent and child GO attributes

%% Read in the new information about parent GO categories

% This grabs strings of text so that each element of B is a word
% E.g. '[Term]', 'id', 'GO:00000054'
B = textread('RawData/goslim_yeast.obo.txt','%s');    
B(156:170)

%%
% Find all the places that a new term starts. A "Term" is a GO attribute
% E.g. id: GO:0000054
TermStarts = find(strcmp(B,'[Term]'));

for i = 1:numel(TermStarts)-1                           % For each Term
    TermInfo = B(TermStarts(i):TermStarts(i+1)-1);      % Grab text regarding this term
    
    Goterm = TermInfo(3);     % The third word in the list is the GO id
                              % E.g. GO:0000054
    
    % We are interested in finding parent information:
    % E.g. GO:0000054 is_a GO:0015931
    %      GO:0000910  has_part GO:0008150
    %
    % All the following phrases in the file indicate a parent
    % relationship between a GO attribute and "parent" attributes
    is_a = find(strcmp(TermInfo,'is_a:'));
    has_part = find(strcmp(TermInfo,'has_part'));
    part_of = find(strcmp(TermInfo,'part_of'));
    occurs_in = find(strcmp(TermInfo,'occurs_in'));
    regulates = find(strcmp(TermInfo,'regulates'));
    negreg = find(strcmp(TermInfo,'negatively_regulates'));
    
    % Collect all incides in TermInf where the parent info shows up
    parents = [is_a;has_part;part_of;occurs_in;regulates;negreg];
    
    % Match the GO attribute we're looking at with its column ID in 
    % the GeneToGO matrix
    childID = find(strcmp(GO_list,char(Goterm)));
    
    % In the descendant-ancestor graph, draw an edge from each child to
    % parent
    for j = 1:numel(parents)
        parent = TermInfo(parents(j)+1);
        parentID = find(strcmp(GO_list,char(parent)));
        gograph(childID,parentID) = 1;
    end
    
end

%% Currently gograph(i,j) = 1 iff attribute j is a parent to attribute i.
% We want gograph(i,j) = 1 iff attribute j is an ancestor of attribute i.

for i = 1:n                                 % for each go attribute,
    AncestorList = find(gograph(i,:));      % get immediate parents
    
    while numel(AncestorList) > 0           
        current = AncestorList(1);          % pop an ancestor off
        AncestorList = AncestorList(2:end); % (so remove it from the Ancestor stack)
        
        gograph(i,current) = 1;             % draw an edge from child to ancestor
        
        newAncestors = find(gograph(current,:));        % add ancestors of ancestor
        AncestorList = union(AncestorList,newAncestors);% to the stack
    end
end
 
%% Use this to get a gene-to-attribute dataset that includes information about attribute 
% categories and subcategories.

[m,n] = size(GeneToGO);

GeneToFullGO = GeneToGO; % start with the base gene-to-attribute info                    

for i = 1:m                                   % For each gene,
    MainGOtags = find(GeneToGO(i,:));         % find the main GO attributes.
    AllAncestors = gograph(MainGOtags,:);     % For those attributes, get the incidator matrix for their ancestors
    [~,J,~] = find(AllAncestors);             % The column of each nonzero entry corresponds to an ancestor of a GO attribute.
    GeneToFullGO(i,unique(J)) = 1;            % So gene i is in some way associated with all these attribute ancestors
    %fprintf('Gene %d has %d more ancestors added. size(AllAncestors) = %d %d \n',i,numel(intersect(unique(J),MainGOtags)),size(AllAncestors,1),size(AllAncestors,2))
end

GeneToFullGO = sparse(GeneToFullGO);

%%
save('MatFiles/AllGOattributes','GeneToFullGO','genes','GO_list');
