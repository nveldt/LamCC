function c = CD_lp_round(Dist)
% Rounds the LP relaxation for cluster deletion into an actual clustering.
% This will return a 4-approximation to the optimal solution, as proven by
% Veldt, Gleich, and Wirth (2017)

n = size(Dist,1);
Vinds = (1:n)';
C = [];

while numel(Vinds >0)
    i = randi(numel(Vinds));          % select random index from Vinds
    pivot = Vinds(i);                 % index for pivot
    scores = full(Dist(pivot,Vinds)); % list of distances from pivot
    
    Tindicator = scores < 1/2;        % it is important in the proof to make this a strict inequality
    T = find(Tindicator);
    
    Tscores = scores(T);
    average = mean(Tscores);
    if average <= 1/4
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
