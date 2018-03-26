% Load plot data and plot the result

load PlotDataFor_Cornell5
name = 'Cornell5';
%% Probability Plot
xs = linspace(.005,.25,20);
Good = True;
Bad = Spur;
clf

lw = 1.5;
ms = 10;
darkgreen = [0 .85 0];
set(gca,'fontsize', 11);

figure(1)
jump = 3;
pts = xs(1:jump:end);
clus = numClus(1:jump:end);
%clus = round(n./clus);
text(pts,-0.040*ones(numel(pts),1),num2str(clus))
hold on

% Student/Faculty
plot(xs, Good(1,:),'.-','linewidth',lw,'markersize',ms,'color','b')
% Dorm
plot(xs, Good(3,:),'.-','linewidth',lw,'markersize',ms,'color','r')
% Graduation Year
plot(xs, Good(4,:),'.-','linewidth',lw,'markersize',ms,'color',darkgreen)

legend('S/F','Dorm','Year')
legend('Location','Best');
legend boxoff

plot(xs, Bad(1,:),'--','linewidth',1,'markersize',10,'color','b')
plot(xs, Bad(3,:),'--','linewidth',1,'markersize',10,'color','r')
plot(xs, Bad(4,:),'--','linewidth',1,'markersize',10,'color',darkgreen)

ylim([-.08,1])

xlabel('\lambda*n','fontsize',16);
ylabel({'Probability of Shared Attribute'});

hold off

%%
set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf(strcat('Plots/final_',name,'.eps')),'-depsc2');
Process_AtendHeader(strcat('Plots/final_',name,'.eps'),'');

%%

set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf(strcat('Plots/nocropfinal_',name,'.eps')),'-depsc2');

