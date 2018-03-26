%% Download and store yeast interection data from the String Database

fileID = fopen('RawData/4932.protein.links.v10.txt','r');

tline = fgetl(fileID);
disp(tline)

% First build a list that goes from indices to ORF names
List = {};  
count = 1;
while ischar(tline)
    count = count + 1;
    tline = fgetl(fileID);
    DotPlaces = strfind(tline,'.');
    Spaces = strfind(tline,' ');   
    gene1 = tline(DotPlaces(1)+1:Spaces(1)-1);
    
    % If we haven't already put it in the list, put it in the list now
    location = find(strcmp(List,gene1));
    if numel(location) == 0
        List = [List; gene1];
    end
   
end

size(List)
fclose(fileID);


%% Now fill the interaction score matrix

fileID = fopen('RawData/4932.protein.links.v10.txt','r');
tline = fgetl(fileID);
n = size(List,1);
A = zeros(n);
while ischar(tline)
    tline = fgetl(fileID);
    DotPlaces = strfind(tline,'.');
    Spaces = strfind(tline,' ');   
    gene1 = tline(DotPlaces(1)+1:Spaces(1)-1);
    gene2 = tline(DotPlaces(2)+1:Spaces(2)-1);
    number = str2num(tline(end-3:end));
    %fprintf('%s %s %d\n',gene1, gene2,number);

    location1 = find(strcmp(List,gene1));
    location2 = find(strcmp(List,gene2));

    A(location1,location2) = number;
    A(location2,location1) = number;

end


fclose(fileID);
size(A)
A(1:10,1:10)

%% Save Data
string_labels = List;
save('MatFiles/StringMatrix','A','string_labels')
