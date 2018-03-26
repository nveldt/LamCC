clear

addpath('../Algs/')
addpath('clusterings/')
% Largest connected component of ca-GrQc
load graphs/cagrqc_conn

% Results from algorithms that aren't LambdaCC-based.
%   Outside software:
%       * Metis: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
%       * Graclus: http://www.cs.utexas.edu/users/dml/Software/graclus.html
%       * Louvain: https://github.com/ayanonagon/mkse312_project/tree/master/Community_BGLL_Matlab
%       * InfoMap: http://www.mapequation.org/)
%   Louvain, recursive densest subraph (which doesn't perform well)
%   and InfoMap
load Nonlcc_connected_ca-GrQc.mat

% Use just two clusters for Graclus
cG = cGrac(:,1) + 1;
Gr = LamRange_objs(A,cG,Lams);

% LambdaCC-based algorithms: GrowCluster and Lambda-Louvain
%   * Slow implementations of both in Algs folder
%   * GenLouvain can be used for Lambda-Louvain: http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
load LccAlgs_connected_ca-GrQc.mat

% Recursive Maximum Clique (RMC).
%   Outside software:
%       * PMC library: http://maximumclique.com/
%       * MaximalCliques: https://github.com/aaronmcdaid/MaximalCliques.git
%       These are called in the julia code (see repo folder Algs/julia)
load c2app_cagrqc_conn
RMC = LamRange_objs(A,c2app,Lams);

% Recursive Maximum Quasi-Clique, computed for density = .6
% 	Outside software:
%       * Quick: https://www.comp.nus.edu.sg/~wongls/projects/pattern-spaces/quick-v1/
load recQclique_conn_caGrQc_pt6.mat
RMQC = LamRange_objs(A,c,Lams);

% The LP bounds for LambdaCC
load cagrqc_LPbounds
Bounds = cagrqc_LPbounds;

% The plot looks similar if you just compare against Lambda-Louvain scores
% Bounds = LLs;

[I,J] = find(A);
EdgeList = [I J];
n = size(A,1);

% Plot style variables
lw = 1.5;
ms = 10;

infocolor = [0 .75 0];
rmqcColor = [.65 0 0];
graColor = [1 .6 0];
graColor = [0 .6 .6];

rmqcColor = [0 .6 .6];
rmqcColor = [.9 .5 0];
graColor = [.6 0 0];

% We plot different values of lambda equidistantly
%   This is a pseudo-log scale: the first half are logarithmic, second
%   half are linear, to highlight behavior for very tiny lambda.
xs = 1:numel(Lams);

figure(1)
plot(xs,Gr'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',graColor)
hold on
plot(xs,Louobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,Infoobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',infocolor)
plot(xs,RMQC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',rmqcColor)
plot(xs,RMC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')

% Not all clustering algorithms "win" in a certain regime
        % GrowCluster does worse than Lambda-Louvain.
        % Recursive densest subgraph doesn't do that well anywhere.
        %   (note that it uses a different definition of "density" than
        %       is used by RMC and RMQC)
% plot(xs,GCs./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 .5 0])
% plot(xs,Denobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','m')

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

% Text for clustering names
fsize = 11;
h=text(6.5,4,'Graclus');
set(h,'Rotation',87,'color',graColor,'fontsize',fsize);
 
h2=text(9.6,3.95,'Louvain');
set(h2,'Rotation',87,'color',[.5 0 .75],'fontsize',fsize);

%h3=text(12.2,2.2,'InfoMap');
h3=text(14,3.8,'InfoMap');
set(h3,'Rotation',84,'color',infocolor,'fontsize',fsize);

%h4=text(16.2,1.7,'RMQC');
h4=text(17.8,2.3,'RMQC');
set(h4,'Rotation',77,'color',rmqcColor,'fontsize',fsize);

h4=text(12.2,1.6,'RMC');
set(h4,'Rotation',-10,'color','b','fontsize',fsize);

h4=text(14,1.1,'LamLouv');
set(h4,'Rotation',0,'color','k','fontsize',fsize);

Places = [1 4 7 10 12 15 18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

set(gca,'YTick',[1 2 3 4])
set(gca,'YTickLabel',{'1' '2' '3' '4'});

xlabel('\lambda','fontsize',15);
ylabel('Ratio to LP bound')
ylim([0.9,5])
xlim([1,19])
box off
hold off

%%
set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/cagrqcLP.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');

%% ARI
load arinmi_cagrqc
figure(2)
plot(xs,ariMat(4,:),'.-','linewidth',lw,'markersize',ms,'color',infocolor)
hold on
plot(xs,ariMat(3,:),'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,ariMat(2,:),'.-','linewidth',lw,'markersize',ms,'color',graColor)
plot(xs,ones(numlam,1),'.-','linewidth',lw,'markersize',ms,'color','k')
plot(xs,ariMat(6,:),'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,ariMat(5,:),'.-','linewidth',lw,'markersize',ms,'color',rmqcColor)

ylim([0,1])
xlim([1,19])
xlabel('\lambda','fontsize',15);
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


%% Plot many runs of graclus

figure(4)
plot(1:numel(Lams),LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')
hold on
for t = 1:2:size(cGrac,2)
    cG = cGrac(:,t) + 1;
    max(cG)
plot(1:numel(Lams),AllObjsGrac(t,1:end)'./Bounds(1:end),'.-','linewidth',lw,'markersize',ms,'color',(1-t/(40))*[1 0 0])
end
xlabel('Lambda');
ylabel('Ratio to LP bound')
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
xlabel('Lambda');
ylabel('Ratio to LP Bound')
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
