marker_colors{1} = [019, 103, 131]/256;
marker_colors{2} = [255,158,002] / 256;
marker_colors{3} = [251,056,071] / 256;
marker_colors{4} = [131,064,028] / 256;
marker_colors{5} = [075,116,178] / 256;
marker_colors{6} = [107,112,092] / 256;

% plotChargeCorr1
% hold on;
plotChargeCorr2
hold on;
plotChargeCorr3

% set(h1, 'MarkerEdgeColor', marker_colors{1});hold on;
set(h2, 'MarkerEdgeColor', marker_colors{3});hold on;
set(h3, 'MarkerEdgeColor', marker_colors{4});hold on;

% set(h1_negative, 'Color', marker_colors{1});
% set(h1_negative, 'MarkerFaceColor',  marker_colors{1});

set(h2_negative, 'Color', marker_colors{3});
set(h2_negative, 'MarkerFaceColor',  marker_colors{3});

set(h3_negative, 'Color', marker_colors{4});
set(h3_negative, 'MarkerFaceColor',  marker_colors{4});

% h_fit1.Color =  marker_colors{6};
h_fit2.Color =  marker_colors{6};
h_fit3.Color =  marker_colors{6};


set(gca,'Children',[h2 h3  h2_negative h3_negative  h_fit2 h_fit3 ])
xlim([0 64/2]);
ylim([1e-6, 1e-2]);

set(gca,'fontsize',20);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$\Delta x$','Interpreter','latex');
ylabel('$D(\Delta x)$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',20);
set(get(gca,'YLabel'),'FontSize',20);
set(get(gca,'XLabel'),'FontName','Arial');
set(get(gca,'YLabel'),'FontName','Arial');


xlim([2 32])
xticks([2,4,8,16,32]);
set(gcf,'position',[1000,1000,400,300]);