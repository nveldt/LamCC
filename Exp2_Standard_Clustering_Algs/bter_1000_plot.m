clear
load c_bter_1000_2app
c2app = c;

addpath('../Algs/')
load bter1000main

% Bounds = LLs;

%% Better Bounds


Bounds = [43.969 94.7035 198.2616 404.24 767.88 1336.1 2185 ...
    2783.14 3077.1 3088.8 2957.900 2635.200 2301.800 1959.400 ...
    1611.500 1259.300 903.000 542.475 180.825]';


RMC = LamRange_objs(A,c2app,Lams);

xs = 1:numel(Lams);
figure(1)
color = 'r';

lw = 1.5;
ms = 10;

infocolor = [0 .75 0];
rmqcColor = [.65 0 0];
graColor = [1 .6 0];
graColor = [0 .6 .6];

rmqcColor = [0 .6 .6];
rmqcColor = [.9 .5 0];
graColor = [.6 0 0];

plot(xs,Gr'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',graColor)
hold on
plot(xs,Louobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',[.5 0 .75])
plot(xs,Infoobj./Bounds,'.-','linewidth',lw,'markersize',ms,'color',infocolor)
plot(xs,RMQC'./Bounds,'.-','linewidth',lw,'markersize',ms,'color',rmqcColor)
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

fsize = 11;
h=text(8,4,'Graclus');
set(h,'Rotation',80,'color',graColor,'fontsize',fsize);

h2=text(12,3.5,'Louvain');
set(h2,'Rotation',75,'color',[.5 0 .75],'fontsize',fsize);

h3=text(15.5,3.5,'InfoMap');
set(h3,'Rotation',75,'color',infocolor,'fontsize',fsize);

h4=text(17.9,3.5,'  RMQC');
set(h4,'Rotation',85,'color',rmqcColor,'fontsize',fsize);

h4=text(8.8,2.1,'RMC');
set(h4,'Rotation',-15,'color','b','fontsize',fsize);

h4=text(15,1.35,'LamLouv');
set(h4,'Rotation',0,'color','k','fontsize',fsize);

xlabel('\lambda','fontsize',15);
ylabel('Ratio to LP Bound')
plotname = strcat('Figures/btrmain.eps');
box off
ylim([0.9,5])
xlim([1,19])

hold off

%% Save the details
set_figure_size([2.25*1.75,1.75*1.75]);
plotname = strcat('Figures/btrmain_new.eps');
print(gcf,sprintf(plotname),'-depsc2');
Process_AtendHeader(plotname,'');

%% ARI

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

