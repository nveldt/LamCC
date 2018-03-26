function [C,c] = gamma_delta_rounding(Dist,gamma,delta)

% GAMMA_DELTA_ROUNDING: this is the LP rounding scheme developed by
% Charikar, Guruswami, and Wirth in their paper Clustering with Qualitative
% Information. 
%
% Outline:
%   1. Pick a random vertex u
%   2. Set T = {i \in V : Dist(u,i) \leq gamma}
%   3. If on average Dist(u,i) is less than delta, form cluster C =
%   T\cup{u}
%   4. Else form cluster C = {u}
%   5. Remove cluster C from graph and repeat from step u until no more
%   vertices need to be clustered
%
% Charikar et al specifically chose Gamma = 1/2 and Delta = 1/4
%
% The algorithm takes as input the values of the already solved LP
% relaxation, Dist

n = size(Dist,1);
Vinds = (1:n)';
C = [];

while numel(Vinds >0)
    i = randi(numel(Vinds));          % select random index from Vinds
    pivot = Vinds(i);                 % index for pivot
    scores = full(Dist(pivot,Vinds)); % list of distances from pivot
    
    Tindicator = scores <= gamma;
    T = find(Tindicator);
    
    Tscores = scores(T);
    average = mean(Tscores);
    if average < delta
        NewCinds = [pivot; Vinds(T)];
    else
        NewCinds = pivot;
    end
    
    ci = zeros(n,1);
    ci(NewCinds) = 1;           % Forms a cluster indicator vector
    C = [C ci];
    Vinds = setdiff(Vinds,NewCinds);    % updates Vinds
end
m = size(C,2);
v = 1:m;
c = C*v';

