%% Load LFR graph, and add path to other algorithms
clear
load A_bter_1000_exp
n = size(A,1);

Kvals = [2:10 15:5:45, 50:10:90, 100:50:700];
addpath('../Community_BGLL_Matlab/')
addpath('Algs/')
addpath('../Lambda_louvain/')
addpath('../find_densest_subgraph-master/matlab_wrapper_loadlibrary/')
addpath('../metis-5.0.2/metismex/')
addpath('/Users/nateveldt/GitHubRepos/LFR_lcc/quick/quick-source')

[I,J] = find(A);
EdgeList = [I,J];
%%
i = 1;
Lams = [logspace(-4,-1,10), .15:.1:.95];
numlam = numel(Lams);
numlams = numlam;
numtimes = 1;
Louobj = zeros(numlams,numtimes);
Metisobj = Louobj;
Graclusobj = Louobj;
Denobj = Louobj;
Infoobj = Louobj;

n = size(A,1);

% Louvain
COMTY = cluster_jl_cpp(triu(full(A)),1,1,0,1);
cLou = COMTY.COM{end};
Louobj(:,i) = LamRange_objs(A,cLou,Lams);

% Infomap
cInfo = Run_InfoMap(A,1);
Infoobj(:,i) = LamRange_objs(A,cInfo,Lams);

% Recursive densest subgraph
cDen = meng_recursive_dense_subgraph(A);
Denobj(:,i) = LamRange_objs(A,cDen,Lams);

% Metis
cMet = zeros(n,numel(Kvals));
for ii = 1:numel(Kvals)
    kk = Kvals(ii);
    c = metismex('PartGraphRecursive',sparse(A),kk);
    cMet(:,ii) = c+1;
end
AllObjsMet = zeros(numel(Kvals),numlam);
for ii = 1:numel(Kvals);
    AllObjsMet(ii,:) = LamRange_objs(A,cMet(:,ii),Lams);
end
[Metisobj(:,i), metisLamclus] = min(AllObjsMet);

% Graclus
cGrac = Run_Graclus(A,EdgeList,Kvals);
AllObjsGrac = zeros(numel(Kvals),numlam);
for ii = 1:numel(Kvals);
    AllObjsGrac(ii,:) = LamRange_objs(A,cGrac(:,ii),Lams);
end
[Graclusobj(:,i), gracLamclus] = min(AllObjsGrac);

%% Recursive maximum quasi-clique
% 
% den = .5;
% minsize = 50;
% cRMQ5 = recursive_maximum_quasiclique(A,den,minsize);

den = .55;
minsize = 50;
cRMQ55 = recursive_maximum_quasiclique(A,den,minsize);

den = .6;
minsize = 50;
cRMQ6 = recursive_maximum_quasiclique(A,den,minsize);
% 
% den = .65;
% minsize = 50;
% cRMQ65 = recursive_maximum_quasiclique(A,den,minsize);
% 
% den = .75;
% minsize = 50;
% cRMQ75 = recursive_maximum_quasiclique(A,den,minsize);
% 
% den = .8;
% minsize = 50;
% cRMQ8 = recursive_maximum_quasiclique(A,den,minsize);
% 
% den = .85;
% minsize = 50;
% cRMQ85 = recursive_maximum_quasiclique(A,den,minsize);
% 
% den = .9;
% minsize = 50;
% cRMQ9 = recursive_maximum_quasiclique(A,den,minsize);
% 
% den = .95;
% minsize = 50;
% cRMQ95 = recursive_maximum_quasiclique(A,den,minsize);

%%
cRMQ = cRMQ6;
RMQC = LamRange_objs(A,cRMQ,Lams);

cG = cGrac(:,1) + 1;
cM = cMet(:,1);

Gr = LamRange_objs(A,cG,Lams);
Me = LamRange_objs(A,cM,Lams);

GrMe = Gr;




%% Grow-Cluster and Lambda-Louvain
addpath('/Users/nateveldt/GitHubRepos/LFR_lcc/Lambda_louvain')
GCclus = zeros(n,numlam);
LLclus = zeros(n,numlam);
GCs = zeros(numlam,1);
LLs = zeros(numlam,1);

for ii = 1:numlam
    cGrow = GrowCluster(A,Lams(ii),1);
    cLL = lambda_louvain(A,Lams(ii),ones(n,1));
    GCs(ii) = lamCCobj(A,Lams(ii),cGrow);
    LLs(ii) = lamCCobj(A,Lams(ii),cLL);
    GCclus(:,ii) = cGrow;
    LLclus(:,ii) = cLL;
end
LLobj = LLs;
GCobj = GCs;

% save('Plotdata/bterLLs')
%save('Plotdata/bter100main')
% Save it
save('/Users/nateveldt/GitHubRepos/LFR_lcc/Code/Plotdata/bter1000main')
%%
Bounds = LLs;

lw = 2;
ms = 15;
xs = 1:numel(Lams);
figure(1)
plot(xs,Infoobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','c')
hold on
plot(xs,Louobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .5])
%plot(xs,Denobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','m')
plot(xs,GrMe'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[1 .6 0])

plot(xs,RMQC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','g')
%plot(xs,GCobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','y')
plot(xs,LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')

Ticks = [];
for ii = 1:2:numel(Lams)
    Ticks = [Ticks; cellstr(num2str(Lams(ii),2))];
end
set(gca,'XTick',1:2:numel(Lams))
set(gca,'XTickLabel',Ticks);
xlabel('Lambda');
ylabel('Ratio to Lambda-Louvain Score')
plotname = strcat('Figures/btrmain.eps');

%title(plotname)
ylim([0.9,5])
xlim([4,19])
legend('InfoMap','Louvain','Graclus','RMQC','Lam-Louv')
legend('Location','North');
legend boxoff
hold off

%% Save the details
cd /Users/nateveldt/GitHubRepos/LFR_lcc/Code
set_figure_size([2.25*1.5,1.75*1.5]);
plotname = strcat('Figures/btrmain.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');