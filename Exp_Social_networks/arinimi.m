function [ari,nmi] = arinimi(info,c)
addpath('../Algs/')
ari = zeros(size(info,2),1);
nmi = zeros(size(info,2),1);

for i = 1:size(info,2)
    ari(i) = ARI(info(:,i)+1,c);
    nmi(i) = NMI(info(:,i)+1,c);
end

end