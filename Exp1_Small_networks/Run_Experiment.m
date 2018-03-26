addpath('../Algs/')
addpath('../Networks/')


for i = 2:6

    if i == 1
        name = 'karate';
    elseif i ==2
        name = 'dolphins';
    elseif i ==3
        name = 'lesmis';
    elseif i ==4
        name = 'polbooks';
    elseif i == 5
        name = 'football';
    else
        name = 'adjnoun';
    end
    
    load(strcat(name,'.mat'))

load(strcat(name,'.mat'));
load(strcat(name,'SparsestCut'));
SC = ScaledSC;
A = spones(Problem.A);

Lams = [logspace(-3,-1,7) .15:.05:.95,.99];

[Bounds, LPs, LPtime,Objectives, Times] = run_all_algs(A,Lams,0);

% Ratios
Piv = Objectives(:,2)./Bounds';
GC = Objectives(:,3)./Bounds';
LP4 = Objectives(:,4)./Bounds';
LP5 = Objectives(:,5)./Bounds';
LP3 = Objectives(:,6)./Bounds';
BG11 = Objectives(:,7)./Bounds';
LL = Objectives(:,8)./Bounds';

Pt = Times(:,2);
GCt = Times(:,3);
LP4t = Times(:,4);
LP5t = Times(:,5);
LP3t = Times(:,6);
BG11t = Times(:,7);
LLt = Times(:,8);

save(strcat('PlotDataFor',name));

end
