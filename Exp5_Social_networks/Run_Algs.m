addpath('~/data/Facebook100/')
addpath('~/GitHubRepos/GenLouvain-2.1/')


Graphs= {'Swarthmore42','Caltech36'};

%%
for i = 1:size(Graphs,2)
    
name = char(Graphs(i));

load(name)


n = size(A,1);
limit = 10;
d = sum(A,2)';
e = ones(1,n);
n = size(A,1);
 
% Student/Fac, Major, Dorm, Graduation year
local_info = local_info(:,[1, 3, 5, 6]);
p = randperm(n);
spurious = local_info(p,:);


% Run degree-weighted LambdaCC

w = d;
Lams = linspace(.005,.25,20)/n;
numlam = numel(Lams);
limit = 1000;
True = [];
Spur = [];
Base = [];
aTrue = [];
nTrue = [];
nSpur = [];
aSpur = [];
numClus = zeros(numlam,1);
next = 1;
for lam = Lams
    %B = @(i) A(:,i) - lam*w'*w(i);
    B = full(A - w'*w*lam);
    c = iterated_genlouvain(B,limit);
    [~, true, base] = ProbwithClustering(A,local_info,c);
    [~, spur,~] = ProbwithClustering(A,spurious,c); 
    
    [aritrue,nmitrue] = arinimi(local_info,c);
    [arispur,nmispur] = arinimi(spurious,c);
    
    fprintf('\nLambda = %f, num com %d, avg com %.2f \n',lam, max(c),(n/max(c)))
    fprintf(' \t \t \t S/F  \t Major \t Dorm \t Grad\n')
    fprintf('Two random share att: %f %f %f %f \n',base(1),base(2), base(3), base(4))
    fprintf('In cluster spurious : %f %f %f %f \n',spur(1),spur(2),spur(3),spur(4))
    fprintf('In cluster share att: %f %f %f %f \n',true(1),true(2),true(3),true(4))
    fprintf('ari spurious : %f %f %f %f \n',arispur(1),arispur(2),arispur(3),arispur(4))
    fprintf('ari true     : %f %f %f %f \n',aritrue(1),aritrue(2),aritrue(3),aritrue(4))
    fprintf('nmi spurious : %f %f %f %f \n',nmispur(1),nmispur(2),nmispur(3),nmispur(4))
    fprintf('nmi true     : %f %f %f %f \n',nmitrue(1),nmitrue(2),nmitrue(3),nmitrue(4))
    numClus(next) = max(c);
    next = next+1;
    True = [True true];
    Spur = [Spur spur];
    Base = [Base base];
    aTrue = [aTrue aritrue];
    nTrue = [nTrue nmitrue];
    aSpur = [aSpur arispur];
    nSpur = [nSpur nmispur];
end

save(strcat('PlotDataFor_',name),'True','Spur','Base','aTrue','nTrue','aSpur','nSpur','Lams','numClus')

end