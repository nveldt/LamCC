%% Graclus plots
% Load a network

clear
addpath('graphfiles/')
addpath('Algs/')
graph = 'scom10k15';    
DisagreeFlag = 0;           % 1 if you want to plot weight of disagreements
DisOrAg = 'agreements';
if DisagreeFlag
    DisOrAg = 'disagreements';
end
MetOrGrac = 'Metis';

load(strcat('graphfiles/',graph))
n = size(A,1);
FigNum = 2;
%% Get clusterings for a range of clusters k


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
    AllObjs(i,:) = LamRange(A,cMany(:,i),Lams,DisagreeFlag);
end

TruthObjs = LamRange(A,truth,Lams,DisagreeFlag);
LouObjs = LamRange(A,cLou,Lams,DisagreeFlag);
InfoMap = LamRange(A,cInfo,Lams,DisagreeFlag);
Giant = LamRange(A,ones(n,1),Lams,DisagreeFlag);
Singletons = LamRange(A,(1:n)',Lams,DisagreeFlag);
Dens = LamRange(A,cDen,Lams,DisagreeFlag);
GCs = zeros(numlam,1);

for i = 1:numlam
    GCs(i) = LamRange(A,FullGClist(:,i),Lams(i),DisagreeFlag);
end

%%
% load(strcat('clusterings/GCexp_',graph))
% GCexps = zeros(numlam,1);
% for i = 1:numlam
%     GCexps(i) = LamRange_objs(A,FullGClist(:,i),Lams(i));
% end

%% Plot them
PlotInformation = strcat('Many_', MetOrGrac, 'plots. Data for graph_',graph,'_avg degree 15._',DisOrAg);
save(strcat('PlotData/',graph,'_many',MetOrGrac,'_',DisOrAg))

%%
low = .875;
figure(FigNum)
loglog(Lams,InfoMap,'linewidth',2,'color','g')
hold on
loglog(Lams,TruthObjs,'o','color','k')
loglog(Lams,LouObjs,'color','y')
loglog(Lams,GCs,'linewidth',2,'color','b')
%loglog(Lams,GCexps,'linewidth',1,'color',[.5 0 .5])
loglog(Lams,Giant,'--','linewidth',2,'color','r')
loglog(Lams,Singletons,'--','linewidth',2,'color','k')
loglog(Lams,Dens,'c')

if DisagreeFlag
    ylim([0,max(max(LouObjs))])
    ylim([0,1e5])
    xlim([1e-5,1])
else
    ylim([low,1.01])
end

for i = 1:2:numel(Kvals)
    loglog(Lams,AllObjs(i,:),'color',i/((numel(Kvals) + 3))*[1 1 1]);
end
