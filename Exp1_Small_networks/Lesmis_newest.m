%% Here we reproduce the plots from the main text of the paper
% Load data run in the Run_experiments.m file
addpath('PlotData/')
load PlotDataForLP3lesmis.mat
L = round(Lams,3);

% To emphasize the behavior of algorithms for very small lambda, we use a pseudo-logarithmic scale
% (first half of axis is logarithmic, second half is plotted linearly)
Lams = 1:numel(Lams);

%% Plot everything for starters
pivcol = [0 .75 0];
gccol = [.6 0 .6];
icmcol = [1 .5 0]*.95;
figure(1); clf;
hax=axes; 
lw = 1.5;
design = '.-';
ms = 12;
plot(Lams,Piv,design,'LineWidth',lw,'markersize',ms,'color',pivcol);
hold on
plot(Lams,BG11,design,'LineWidth',lw,'markersize',ms,'color',icmcol);
plot(Lams,LL,design,'LineWidth',lw,'markersize',ms,'color','k');
plot(Lams,GC,design,'LineWidth',lw,'markersize',ms,'color',gccol);
plot(Lams,LP3,design,'LineWidth',lw,'markersize',ms,'color','b');

% Text for stuff
set(gca,'fontsize', 11);

h=text(12.8,1.01,'ICM');
fsize = 12;
set(h,'Rotation',0,'color',icmcol,'fontsize',fsize);

h=text(6.8,1.9,'Pivot');
set(h,'Rotation',0,'color',pivcol,'fontsize',fsize);

h=text(17,1.95,'threeLP');
set(h,'Rotation',0,'color','b','fontsize',fsize);

h=text(8.5,1.55,'GrowClus');
set(h,'Rotation',0,'color',gccol,'fontsize',fsize);

h=text(19,1.08,'LamLouv');
set(h,'Rotation',0,'color','k','fontsize',fsize);

set(gca,'YTick',[1 1.5 2 2.5])
set(gca,'YTickLabel',{'1' '1.5' '2' '2.5'});


%% Find out where the sparsest cut line goes
darkyellow = [.9 .9 0]
lamcolor = [0 .6 .6]
lamcolor2 = [.0 .6 .6]*.75
for j = 1:numel(L)
    if SC > L(j)
        here = j;
    end
end
topy = 2.35;
place = (2*here+1)/2;
plot([place place],[0 topy],'--','LineWidth',lw,'color',lamcolor)

h=text(place,.68,'\lambda*');
set(h,'Rotation',0,'color',lamcolor2,'fontsize',fsize+1);


box off
set(gca,'XTick',[4,6,10,14,18,22])

set(gca,'XTickLabel',{'.01','.05','.25','.5','.7','.9'});
xlim([2,numel(Lams)]);

ylim([.9,topy]);
% legend('Piv','ICM','Lam-Louv','GC','threeLP','\lambda^*')
% legend('Location','North');
% legend boxoff
xlabel('\lambda','fontsize',15);
ylabel({'Ratio to LP Bound'});
hold off
%% Save plot if desired
addpath('../Algs')
set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf(strcat('Figures/newestplot_',name,'.eps')),'-depsc2');
Process_AtendHeader(strcat('Figures/newestplot_',name,'.eps'),'');