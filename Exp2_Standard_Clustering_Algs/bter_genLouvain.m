clear
load c_bter_1000_2app
c2app = c;

addpath('Algs/')
load Plotdata/bter1000main

Bounds = LLs;

%% Compute the lambda-louvain clustering a different way
addpath('../../GenLouvain-2.1/')
LLclus = zeros(n,numlam);
LLs = zeros(numlam,1);
e = ones(1,n);
limit = 3;
for ii = 1:numlam
    %cLL = lambda_louvain(A,Lams(ii),ones(n,1),5);
    lam = Lams(ii);
    B = @(i) A(:,i) - lam*e'*e(i);
    cLL = iterated_genlouvain(B,limit);
    LLs(ii) = lamCCobj(A,Lams(ii),cLL);
    LLclus(:,ii) = cLL;
end
LLobj(:,i) = LLs;


%% Better Bounds


Bounds = [43.969 94.7035 198.2616 404.24 767.88 1336.1 2185 ...
    2783.14 3077.1 3088.8 2957.900 2635.200 2301.800 1959.400 ...
    1611.500 1259.300 903.000 542.475 180.825]';


RMC = LamRange_objs(A,c2app,Lams);

xs = 1:numel(Lams);
figure(1)

lw = 2;
ms = 12;
plot(xs,GrMe'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[1 .5 0])
hold on
%plot(xs,Denobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','m')
plot(xs,Louobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,Infoobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','g')
plot(xs,RMQC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','c')
plot(xs,RMC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')


Ticks = [];
for ii = 1:2:numel(Lams)
    Ticks = [Ticks; cellstr(num2str(Lams(ii),2))];
end
set(gca,'XTick',1:2:numel(Lams))
set(gca,'XTickLabel',Ticks);

Ticks = {'1e-4', '1e-3', '1e-2','.1', '0.25','0.55','0.85'};
Places = [1,4,7,10,12,15,18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

h=text(8,4,'Graclus');
set(h,'Rotation',80,'color',[1 .6 0],'fontsize',9);

h2=text(12,3.5,'Louvain');
set(h2,'Rotation',75,'color',[.5 0 .75],'fontsize',9);

h3=text(15.5,3.5,'InfoMap');
set(h3,'Rotation',75,'color','g','fontsize',9);

h4=text(17.9,3.5,'  RMQC');
set(h4,'Rotation',85,'color','c','fontsize',9);

h4=text(8.8,2.1,'RMC');
set(h4,'Rotation',-15,'color','b','fontsize',9);

h4=text(15,1.35,'Lam-Louv');
set(h4,'Rotation',0,'color','k','fontsize',9);

xlabel('Lambda');
ylabel('Ratio to LP Bound')
plotname = strcat('Figures/btrmain.eps');
box off
%title(plotname)
ylim([0.9,5])
xlim([1,19])
% legend('Graclus','Louvain','InfoMap','RMQC','RMC','Lam-Louv')
% legend('Location','North');
% legend boxoff
hold off

%% Save the details
cd /Users/nateveldt/GitHubRepos/LFR_lcc/Code
set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/btrmain.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');


%% NMI and ARI computations
numalgs = 6;
numlam = numel(Lams);

ariMat = zeros(numalgs,numlam);
nmiMat = zeros(numalgs,numlam);

Clusterings = [cM cG cLou' cInfo cRMQ6 c2app];

for i = 1:numlam
   
    for j = 1:numalgs
    ariMat(j,i) = ARI(Clusterings(:,j),LLclus(:,i));
    nmiMat(j,i) = NMI(Clusterings(:,j),LLclus(:,i));
    end
end


%% ARI
figure(2)
plot(xs,ariMat(4,:),'.-','linewidth',lw,'markersize',ms,'color','g')
hold on
plot(xs,ariMat(3,:),'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,ariMat(2,:),'.-','linewidth',lw,'markersize',ms,'color',[1 .6 0])
%plot(xs,ariMat(1,:),'.-','linewidth',lw,'markersize',ms,'color','y')
plot(xs,ones(numlam,1),'.-','linewidth',lw,'markersize',ms,'color','k')
plot(xs,ariMat(6,:),'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,ariMat(5,:),'.-','linewidth',lw,'markersize',ms,'color','c')

ylim([0,1])
xlim([1,19])
xlabel('Lambda');
ylabel('ARI score with Lam-Louv')

Ticks = {'1e-4', '1e-3', '1e-2','.1', '0.25','0.55','0.85'};
Places = [1,4,7,10,12,15,18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);



hold off
box off

set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/bter_ari.eps');
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

Ticks = {'1e-4', '1e-3', '1e-2','.1', '0.25','0.55','0.85'};
Places = [1,4,7,10,12,15,18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

box off
hold off
set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/bter_nmi.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');

%% Ratio to LL

figure(4)

lw = 2;
ms = 12;
Bounds = LLs;
plot(xs,GrMe'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[1 .5 0])
hold on
%plot(xs,Denobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','m')
plot(xs,Louobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,Infoobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color','g')
plot(xs,RMQC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','c')
plot(xs,RMC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color','b')
plot(xs,LLs./Bounds,'.-','linewidth',lw,'markersize',ms,'color','k')


Ticks = [];
for ii = 1:2:numel(Lams)
    Ticks = [Ticks; cellstr(num2str(Lams(ii),2))];
end
set(gca,'XTick',1:2:numel(Lams))
set(gca,'XTickLabel',Ticks);

Ticks = {'1e-4', '1e-3', '1e-2','.1', '0.25','0.55','0.85'};
Places = [1,4,7,10,12,15,18];
set(gca,'XTick',Places)
set(gca,'XTickLabel',Ticks);

h=text(8,3.7,'Graclus');
set(h,'Rotation',80,'color',[1 .6 0],'fontsize',9);

h2=text(13,3.5,'Louvain');
set(h2,'Rotation',75,'color',[.5 0 .75],'fontsize',9);

h3=text(16.5,3,'InfoMap');
set(h3,'Rotation',80,'color','g','fontsize',9);

h4=text(17.7,2,'  RMQC');
set(h4,'Rotation',78,'color','c','fontsize',9);

h4=text(9.4,1.8,'RMC');
set(h4,'Rotation',-15,'color','b','fontsize',9);

h4=text(15,.85,'Lam-Louv');
set(h4,'Rotation',0,'color','k','fontsize',9);

xlabel('Lambda');
ylabel('Ratio to Lam-Louv')
plotname = strcat('Figures/btrmain.eps');
box off
ylim([0.7,5])
xlim([1,19])
hold off

set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/bter_ll.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');