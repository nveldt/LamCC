% Load plot data and plot the result

load PlotDataFor_Swarthmore42
name = 'Swarthmore42';
%% Probability Plot
xs = linspace(.005,.25,20);
Good = True;
Bad = Spur;
clf
figure(1)
jump = 3;
pts = xs(1:jump:end);
clus = numClus(1:jump:end);
%clus = round(n./clus);
text(pts,0.015*ones(numel(pts),1),num2str(clus))
hold on

% Student/Faculty
plot(xs, Good(1,:),'.-','linewidth',2,'markersize',10,'color','b')
% Dorm
plot(xs, Good(3,:),'.-','linewidth',2,'markersize',10,'color','r')
% Graduation Year
plot(xs, Good(4,:),'.-','linewidth',2,'markersize',10,'color','g')

legend('S/F','Dorm','Year')
legend('Location','Best');
legend boxoff

plot(xs, Bad(1,:),'--','linewidth',1,'markersize',10,'color','b')
plot(xs, Bad(3,:),'--','linewidth',1,'markersize',10,'color','r')
plot(xs, Bad(4,:),'--','linewidth',1,'markersize',10,'color','g')

ylim([-.02,1])

xlabel('n*Lambda');
ylabel({'Probability of Shared Attribute'});

hold off


%%
set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf(strcat('Plots/dLcc_Probs_',name,'.eps')),'-depsc2');
Process_AtendHeader(strcat('Plots/dLcc_Probs_',name,'.eps'),'');

