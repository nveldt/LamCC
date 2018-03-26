%% Graclus plots
% Load a network

clear
addpath('graphfiles/')
addpath('Algs/')
graph = 'scom10k15';    
DisagreeFlag = 1;           % 1 if you want to plot weight of disagreements
load(strcat('graphfiles/',graph))
n = size(A,1);
%% Get clusterings for a range of clusters k

MetOrGrac = 'Metis';
Kvals = [2 3 5 10 100 300 750 1000 2000 3000];
if strcmp(MetOrGrac,'Metis');
    cMany = zeros(n,numel(Kvals));
    for i = 1:numel(Kvals)
        k = Kvals(i);
        c = metismex('PartGraphRecursive',sparse(A),k);
        cMany(:,i) = c+1;
    end
else
cMany = Run_Graclus(A,EdgeList,Kvals);
end
%% Load other clusterings

load(strcat('clusterings/GCs_',graph))
load(strcat('clusterings/cLou_',graph))
load(strcat('clusterings/cInfo_',graph))
load(strcat('clusterings/cDen_',graph))
%% Get objective scores for them

numlam = numel(Lams);
AllObjs = zeros(numel(Kvals),numlam);

for i = 1:numel(Kvals);
    AllObjs(i,:) = LamRange_objs(A,cMany(:,i),Lams);
end

TruthObjs = LamRange_objs(A,truth,Lams);
LouObjs = LamRange_objs(A,cLou,Lams);
InfoMap = LamRange_objs(A,cInfo,Lams);
Giant = LamRange_objs(A,ones(n,1),Lams);
Singletons = LamRange_objs(A,(1:n)',Lams);
Dens = LamRange_objs(A,cDen,Lams);

GCs = zeros(numlam,1);

for i = 1:numlam
    GCs(i) = LamRange_objs(A,FullGClist(:,i),Lams(i));
end

%%

load(strcat('clusterings/GCs_',graph))
GCexps = zeros(numlam,1);
for i = 1:numlam
    GCexps(i) = LamRange_objs(A,FullGClist(:,i),Lams(i));
end
%% Plot them
PlotInformation = strcat('Many ', MetOrGrac, 'plots. Data for 10k node graph with lots of small communities, avg degree 15. Disagreements');
save(strcat('PlotData/scom10k15_many',MetOrGrac,'mistakes'))

%%
figure(1)
title('Graclus Partitionings for a range of k')
loglog(Lams,InfoMap,'linewidth',2,'color','g')
%loglog(Lams,TruthObjs,'o','color','k')
hold on
%loglog(Lams,LouObjs,'color','g')
loglog(Lams,GCs,'linewidth',2,'color','b')
loglog(Lams,Giant,'--','linewidth',2,'color','r')
loglog(Lams,Singletons,'--','linewidth',2,'color','k')
loglog(Lams,Dens,'c')
ylim([0,max(max(LouObjs))])

ylim([0,1e5])
for i = 1:numel(Kvals)
    loglog(Lams,AllObjs(i,:),'color',i/((numel(Kvals) + 3))*[1 1 1]);
end


%% Now plot the weight of agreements divided by total possible agreements

numlam = numel(Lams);
AllObjs = zeros(numel(Kvals),numlam);

for i = 1:numel(Kvals);
    AllObjs(i,:) = LamRange_agreements(A,cMany(:,i),Lams);
end
%%
TruthObjs = LamRange_agreements(A,truth,Lams);
LouObjs = LamRange_agreements(A,cLou,Lams);
InfoMap = LamRange_agreements(A,cInfo,Lams);
Giant = LamRange_agreements(A,ones(n,1),Lams);
Singletons = LamRange_agreements(A,(1:n)',Lams);
GCs = zeros(numlam,1);

for i = 1:numlam
    GCs(i) = LamRange_agreements(A,FullGClist(:,i),Lams(i));
end

%% Plot them
PlotInformation = 'plot data for 10k node graph with lots of small communities, avg degree 15, lots of Graclus clusterings, MAXAGREE plots';
save('PlotData/scom10k15_manyGraclus_agree')

%%
low = .85;
figure(1)
loglog(Lams,InfoMap,'linewidth',2,'color','g')
%loglog(Lams,TruthObjs,'o','color','k')
hold on
%loglog(Lams,LouObjs,'color','g')
loglog(Lams,GCs,'linewidth',2,'color','b')
loglog(Lams,Giant,'--','linewidth',2,'color','r')
loglog(Lams,Singletons,'--','linewidth',2,'color','k')
loglog(Lams,Dens,'c')
ylim([low,1.01])
for i = 1:numel(Kvals)
    loglog(Lams,AllObjs(i,:),'color',i/((numel(Kvals) + 3))*[1 1 1]);
end
title('Graclus Partitionings for a range of k')