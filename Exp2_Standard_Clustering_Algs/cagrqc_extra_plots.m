% NMI, ARI, many metis, and many graclus
clear

cd('/Users/nateveldt/GitHubRepos/LFR_lcc/Code')
addpath('Algs/')
load cagrqc_conn
load clusterings/Nonlcc_connected_ca-GrQc.mat
load clusterings/LccAlgs_connected_ca-GrQc.mat

load c2app_cagrqc_conn
RMC = LamRange_objs(A,c2app,Lams);

load clusterings/recQclique_conn_caGrQc_pt6.mat
cqc = c;
RMQC = LamRange_objs(A,c,Lams);
xlow = 1;
n = size(A,1);

% Use LLs as bounds when they are not available
Bounds = LLs;

% I'm going to get just 2 clusters from cGrac
cG = cGrac(:,1) + 1;
cM = cMet(:,1);
Gr = LamRange_objs(A,cG,Lams);
Me = LamRange_objs(A,cM,Lams);
GrMe = Gr;

lw = 2;
ms = 12;
xs = 1:numel(Lams);

range = [1:19];
Bounds = Bounds(range);
xs = 1:numel(Lams);
figure(1)
plot(xs,Gr(range)'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[1 .6 0])
hold on
%plot(xs,Me(range)'./Bounds,'.-','linewidth',2,'markersize',ms,'color','y')
plot(xs,Louobj(range)./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,Infoobj(range)./Bounds,'.-','linewidth',lw,'markersize',ms,'color','g')
%plot(xs,GCobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 .5 0])
plot(xs,RMQC(range)'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','c')
plot(xs,RMC(range)'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,LLs(range)./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')

Ticks = {'1e-05','0.00022','0.0046','0.1','0.25','0.55','0.85'};

% Text for stuff
h=text(6.5,4,'Graclus');
set(h,'Rotation',90,'color',[1 .6 0],'fontsize',9);

% h=text(5.5,2,'Metis');
% set(h,'Rotation',80,'color','y','fontsize',9);
 
h2=text(9.5,3.8,'Louvain');
set(h2,'Rotation',90,'color',[.5 0 .75],'fontsize',9);

h3=text(12.5,2,'InfoMap');
set(h3,'Rotation',75,'color','g','fontsize',9);

h4=text(17.1,1.6,'  RMQC');
set(h4,'Rotation',75,'color','c','fontsize',9);

h4=text(12.1,1.4,'RMC');
set(h4,'Rotation',-10,'color','b','fontsize',9);

h4=text(14,.85,'Lam-Louv');
set(h4,'Rotation',0,'color','k','fontsize',9);

Places = [1 4 7 10 12 15 18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

set(gca,'YTick',[1 2 3 4])
set(gca,'YTickLabel',{'1' '2' '3' '4'});

xlabel('Lambda');
ylabel('Ratio to Lam-Louv')
ylim([0.7,5])
xlim([xlow,19])
box off
hold off

%
set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/cagrqc_extra.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');

%% NMI and ARI computations
% numalgs = 6;
% numlam = numel(Lams);
% 
% ariMat = zeros(numalgs,numlam);
% nmiMat = zeros(numalgs,numlam);
% 
% Clusterings = [cM cG cLou' cInfo cqc c2app];
% 
% for i = 1:numlam
%    
%     for j = 1:numalgs
%     ariMat(j,i) = ARI(Clusterings(:,j),LLclus(:,i));
%     nmiMat(j,i) = NMI(Clusterings(:,j),LLclus(:,i));
%     end
% end


%% ARI
load arinmi_cagrqc
figure(2)
plot(xs,ariMat(4,:),'.-','linewidth',lw,'markersize',ms,'color','g')
hold on
plot(xs,ariMat(3,:),'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,ariMat(2,:),'.-','linewidth',lw,'markersize',ms,'color',[1 .6 0])
%plot(xs,ariMat(1,:),'.-','linewidth',lw,'markersize',ms,'color','y')
plot(xs,ones(numlam,1),'.-','linewidth',lw,'markersize',ms,'color','k')
plot(xs,ariMat(6,:),'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,ariMat(5,:),'.-','linewidth',lw,'markersize',ms,'color','r')

ylim([0,1])
xlim([1,19])
xlabel('Lambda');
ylabel('ARI score with Lam-Louv')

Ticks = {'1e-05','0.00022','0.0046','0.1','0.25','0.55','0.85'};

Places = [1 4 7 10 12 15 18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

% set(gca,'YTick',[1 2 3 4])
% set(gca,'YTickLabel',{'1' '2' '3' '4'});

hold off
box off

set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/cagrqc_ari.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');

%% And now NMI
figure(3)
plot(xs,nmiMat(4,:),'.-','linewidth',lw,'markersize',ms,'color','g')
hold on
plot(xs,nmiMat(3,:),'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,nmiMat(1,:),'.-','linewidth',lw,'markersize',ms,'color',[1 .6 0])
plot(xs,ones(numlam,1),'.-','linewidth',lw,'markersize',ms,'color','k')
plot(xs,nmiMat(6,:),'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,nmiMat(5,:),'.-','linewidth',lw,'markersize',ms,'color','c')
%ylim([0.9,5])
xlim([1,19])
xlabel('Lambda');
ylabel('NMI score with Lam-Louv')
Ticks = {'1e-05','0.00022','0.0046','0.1','0.25','0.55','0.85'};

Places = [1 4 7 10 12 15 18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

% set(gca,'YTick',[1 2 3 4])
% set(gca,'YTickLabel',{'1' '2' '3' '4'});
box off
hold off
set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/cagrqc_nmi.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');

%% Plot many runs of graclus

figure(4)
plot(1:numel(Lams),LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')
hold on
for t = 1:2:size(cGrac,2)
    cG = cGrac(:,t) + 1;
    max(cG)
plot(1:numel(Lams),AllObjsGrac(t,1:end)'./Bounds(1:end),'.-','linewidth',lw,'markersize',ms,'color',(1-t/(40))*[1 .5 0])
end
xlabel('\lambda','fontsize',15);
ylabel('Ratio to Lam-Louv')
Places = [1 4 7 10 12 15 18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

set(gca,'YTick',[1 2 3 4])
set(gca,'YTickLabel',{'1' '2' '3' '4'});
ylim([0.8,5])
xlim([1,19])
box off
hold off

set_figure_size([2.25*2,1.75*2]);
plotname = strcat('Figures/cagrqc_manygraclus.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');


%% Plot many runs of densest subgraph, to show that the quasi-clique clusterings converge to the 2app for CD
rmqc = zeros(10,numlams);

load clusterings/recQclique_conn_caGrQc_pt5.mat
rmqc(3,:) = LamRange_objs(A,c,Lams);

load clusterings/recQclique_conn_caGrQc_pt55.mat
rmqc(4,:) = LamRange_objs(A,c,Lams);

load clusterings/recQclique_conn_caGrQc_pt6.mat
rmqc(5,:) = LamRange_objs(A,c,Lams);
load clusterings/recQclique_conn_caGrQc_pt65.mat
rmqc(6,:) = LamRange_objs(A,c,Lams);

load clusterings/recQclique_conn_caGrQc_pt7.mat
rmqc(7,:) = LamRange_objs(A,c,Lams);

load clusterings/recQclique_conn_caGrQc_pt75.mat
rmqc(8,:) = LamRange_objs(A,c,Lams);

load clusterings/recQclique_conn_caGrQc_pt8.mat
rmqc(9,:) = LamRange_objs(A,c,Lams);


load clusterings/recQclique_conn_caGrQc_pt85.mat
rmqc(10,:) = LamRange_objs(A,c,Lams);



figure(5)
plot(1:numel(Lams),LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')
hold on
for t = 1:10
plot(1:numel(Lams),rmqc(t,:)'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',(1-(t/(10)))/1.3*[0 1 1])
hold on
end
plot(1:numel(Lams),RMC'./Bounds,'.-','linewidth',2,'markersize',ms,'color','b')
hold off
xlabel('\lambda','fontsize',15);
ylabel('Ratio to Lam-Louv')
Places = [1 4 7 10 12 15 18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

set(gca,'YTick',[1 2 3 4])
set(gca,'YTickLabel',{'1' '2' '3' '4'});
ylim([0.8,5])
xlim([1,19])
box off
hold off

set_figure_size([2.25*2,1.75*2]);
plotname = strcat('Figures/cagrqc_manydense.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');
