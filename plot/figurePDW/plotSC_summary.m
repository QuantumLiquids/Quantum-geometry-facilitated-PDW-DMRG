
plotSC1
hold on;
plotSC2
plotSC3
marker_colors{1} = [019, 103, 131]/256;
marker_colors{2} = [255,158,002] / 256;
marker_colors{3} = [251,056,071] / 256;
marker_colors{4} = [131,064,028] / 256;
marker_colors{5} = [075,116,178] / 256;
marker_colors{6} = [107,112,092] / 256;


ylim([1e-3, 1e-1]);
xlim([2 32])
xticks([2,4,8,16,32]);

set(h1, 'MarkerEdgeColor', marker_colors{1});hold on;
set(h2, 'MarkerEdgeColor', marker_colors{3});hold on;
set(h3, 'MarkerEdgeColor', marker_colors{4});hold on;

h_fit1.Color =  marker_colors{6};
h_fit2.Color =  marker_colors{6};
h_fit3.Color =  marker_colors{6};

set(gca,'Children',[h1 h2 h3 h_fit1 h_fit2 h_fit3 ])
% Set the labels and title
set(gca,'fontsize',20);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$\Delta x$','Interpreter','latex');
%Phi(x) = \langle\Delta(0)^\dagger \Delta(x)\rangle
ylabel('$\Phi(\Delta x) \cdot (-1)^{\Delta x}$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',20);
set(get(gca,'YLabel'),'FontSize',20);
set(get(gca,'XLabel'),'FontName','Arial');
set(get(gca,'YLabel'),'FontName','Arial');

legend_entries = {'3-leg \delta=1/8', '4-leg \delta=1/32', '4-leg \delta=1/8'};
l=legend([h1,h2,h3],legend_entries, 'Location', 'best');
set(l,'Box','off');
% set(l,'Interpreter','latex');
set(l,'Fontsize',20);
set(l,'Location','SouthWest');

set(gcf,'position',[1000,1000,400,300]);