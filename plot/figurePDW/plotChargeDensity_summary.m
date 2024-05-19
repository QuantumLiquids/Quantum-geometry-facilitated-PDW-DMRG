marker_colors{1} = [019, 103, 131]/256;
marker_colors{2} = [255,158,002] / 256;
marker_colors{3} = [251,056,071] / 256;
marker_colors{4} = [131,064,028] / 256;
marker_colors{5} = [075,116,178] / 256;
marker_colors{6} = [107,112,092] / 256;

plotChargeDensity1
hold on;
plotChargeDensity2
plotChargeDensity3

set(hA1, 'MarkerEdgeColor', marker_colors{1});hold on;
set(hA2, 'MarkerEdgeColor', marker_colors{1});hold on;
set(hB, 'MarkerEdgeColor', marker_colors{3});hold on;
set(hC, 'MarkerEdgeColor', marker_colors{4});hold on;


h_fit_A1.Color =  marker_colors{6};
h_fit_A2.Color =  marker_colors{6};
h_fit_B.Color =  marker_colors{6};
h_fit_C.Color =  marker_colors{6};

set(gca,'Children',[hA1 hA2 hB hC h_fit_A1 h_fit_A2 h_fit_B h_fit_C])
xlim([0 64/2]);
% ylim([1e-3, 1e-1]);
% xlim([2 32])
% xticks([2,4,8,16,32]);

set(gca,'fontsize',20);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$ x$','Interpreter','latex');
ylabel('$\langle n(x)\rangle$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',20);
set(get(gca,'YLabel'),'FontSize',20);
set(get(gca,'XLabel'),'FontName','Arial');
set(get(gca,'YLabel'),'FontName','Arial');



% l = legend([hA1, hA2, hB, hC], {'Parameter A, s-orbital', 'Parameter A, d-orbital', 'Parameter B', 'C'}, 'Location', 'best');
% 
% set(l,'Box','off');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');

% legend_entries = {'Parameter A', 'B', 'C'};
% l=legend([h1,h2,h3],legend_entries, 'Location', 'best');
% set(l,'Box','off');
% % set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');