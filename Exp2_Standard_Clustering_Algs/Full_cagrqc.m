clear

cd('/Users/nateveldt/GitHubRepos/LFR_lcc/Code')
addpath('Algs/')
load clusterings/Nonlcc_connected_ca-GrQc.mat
load clusterings/LccAlgs_connected_ca-GrQc.mat
load ~/data/MatlabGraphs/ca-GrQc.mat

load cagrqc_conn

size(A)
[I,J] = find(A);
EdgeList = [I J];
n = size(A,1);
color = 'r';
% Use LLs as bounds when they are not available
Bounds = LLs;

% still unsure of first (.002) and third (.0016)
Bounds = [LLs(1),233.7,617.5,1114, 1396,1610,1954.62,2444.75,3064.449,3580.458,3598.123,3358.94,2999.05,2585.05,2154.85,1701.9,1231.63,742.95,250.00]';

% Cluster deletion 2-approximation
load c2app_cagrqc_conn
RMC = LamRange_objs(A,c2app,Lams);

% Get the density clustering
load clusterings/recQclique_conn_caGrQc_pt6.mat
RMQC = LamRange_objs(A,c,Lams);

% I'm going to get just 2 clusters from cGrac
cG = cGrac(:,1) + 1;
cM = cMet(:,1);
Gr = LamRange_objs(A,cG,Lams);
Me = LamRange_objs(A,cM,Lams);
GrMe = Gr;

lw = 2;
ms = 12;
xs = 1:numel(Lams);

figure(1)
plot(xs,GrMe'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[1 .6 0])
hold on
plot(xs,Louobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,Infoobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','g')
%plot(xs,GCobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 .5 0])
plot(xs,RMQC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',color)
plot(xs,RMC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','b')
%plot(xs,LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')

Ticks = [];

spacing = 3;
for ii = 1:spacing:numel(Lams)
    Ticks = [Ticks; cellstr(num2str(Lams(ii),2))];
end
set(gca,'XTick',1:spacing:numel(Lams))
set(gca,'XTickLabel',Ticks);


Ticks = {'1e-05'
    '0.00022'
    '0.0046'
    '0.1'
    '0.25'
    '0.55'
    '0.85'};

% Text for stuff
h=text(6.4,4,'Graclus');
set(h,'Rotation',90,'color',[1 .6 0],'fontsize',9);
 
h2=text(9.4,3.8,'Louvain');
set(h2,'Rotation',90,'color',[.5 0 .75],'fontsize',9);

h3=text(12,2,'InfoMap');
set(h3,'Rotation',75,'color','g','fontsize',9);

h4=text(15.9,1.6,'  RMQC');
set(h4,'Rotation',35,'color',color,'fontsize',9);

h4=text(12.2,1.6,'RMC');
set(h4,'Rotation',-10,'color','b','fontsize',9);

% h4=text(14,1.1,'Lam-Louv');
% set(h4,'Rotation',0,'color','k','fontsize',9);

Places = [1 4 7 10 12 15 18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

set(gca,'YTick',[1 2 3 4])
set(gca,'YTickLabel',{'1' '2' '3' '4'});

xlabel('Lambda');
ylabel('Ratio to LP bound')
ylim([0.9,5])
xlim([1,19])
box off
% legend('Graclus','Louvain','InfoMap','RMQC','RMC','Lam-Louv')
% legend('Location','Best')
% legend boxoff
hold off


set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/caGrQc_fullLPresults.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');


