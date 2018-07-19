function objs = LamRange_objs(A,c,Lams)
% Given a network and a clustering on it, and a range of Lambda values,
% calculate the weight of mistakes for each value of Lambda. Output a list
% of objective values in a vector the same length as Lams

n = size(A,1);
assert(numel(c) == n);

if min(c) == 0
    c = c+1;
end

[pos, neg] = Mistakes(A,c);

objs = Lams;
for i = 1:numel(Lams);
    lam = Lams(i);
    objs(i) = (1-lam)*pos + lam*neg;
end