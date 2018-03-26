% Read a .txt file whose first line is the name of a network
% Second line is the number of nodes
% Third line is the algorithm that was run
% Fourth line is algorithm details or inputs
% Fifth is the time it took the algorithm to run
% The rest is an integer per line telling you which cluster each node
% belongs to.
%
% This is exactly the format we use to store a clustering generating by
% running our Julia code which repeatedly calls the PMC library

function Clustering_txt2mat(filename)

  file = fopen(filename,'r');

Network = fgetl(file);
n = fgetl(file);
alg = fgetl(file);
details = fgetl(file);
runtime = fgetl(file);

c = fscanf(file,'%d');

save(strcat(filename,'.mat'),'c','alg','details','n','Network','runtime')