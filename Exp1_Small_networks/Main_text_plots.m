%% Here we reproduce the plots from the main text of the paper
% Load data run in the Run_experiments.m file
addpath('PlotData/')
load PlotDataForLP3football.mat
L = round(Lams,3);

% To emphasize the behavior of algorithms for very small lambda, we use a pseudo-logarithmic scale
% (first half of axis is logarithmic, second half is plotted linearly)
Lams = 1:numel(Lams);

%% Plot everything for starters

color = [1 .6 0];
figure(1); clf;
hax=axes; 
lw = 2;
design = '.-';
ms = 15;
plot(Lams,Piv,design,'LineWidth',lw,'markersize',ms,'color','g');
hold on
plot(Lams,BG11,design,'LineWidth',lw,'markersize',ms,'color',color);
plot(Lams,LL,design,'LineWidth',lw,'markersize',ms,'color','k');
plot(Lams,GC,design,'LineWidth',lw,'markersize',ms,'color','r');
plot(Lams,LP3,design,'LineWidth',lw,'markersize',ms,'color','b');

% In practice LP5 and LP3 tend to look almost identical
%plot(Lams,LP5,'.-','LineWidth',2,'markersize',15,'color','m');

%% Find out where the sparsest cut line goes

for j = 1:numel(L)
    if SC > L(j)
        here = j;
    end
end
topy = 2.95;
place = (2*here+1)/2;
plot([place place],[0 topy],'LineWidth',lw,'color','y')
box off
set(gca,'XTick',[4,6,10,14,18,22])

set(gca,'XTickLabel',{'.01','.05','.25','.5','.7','.9'});
xlim([2,numel(Lams)]);

% Uncomment this if you want linear axes
%xlim([0,1]);

ylim([.9,topy]);
legend('Piv','ICM','Lam-Louv','GC','threeLP','\lambda^*')
legend('Location','North');
legend boxoff
xlabel('Lambda');
ylabel({'Ratio to LP bound'});
hold off
%% Save plot if desired
set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf(strcat('Figures/plot_',name,'.eps')),'-depsc2');
Process_AtendHeader(strcat('Figures/plot_',name,'.eps'),'');