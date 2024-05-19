U = 8;
colororder("gem");
C = colororder;

% Define marker size
my_marker_size = 100;

% fully gapped
Usd = [2,4,6,8];
delta = [0,0,0,0];
h1 = scatter(Usd, delta, my_marker_size, "filled"); hold on;

Usd = [8];
delta = [0];
scatter(Usd, delta, my_marker_size, "filled", 's', 'MarkerFaceColor', h1.CData); hold on;

% (pi,pi) CDW
Usd = [0,1];
delta = [0,1/32];
h2 = scatter(Usd, delta, my_marker_size, "filled"); hold on;

Usd = [0.3, 0.3];
delta = [0,1- 0.9863403809424479];
% Use the same color as h2 for h2s and set marker size to my_marker_size
h2s = scatter(Usd, delta, my_marker_size, "filled", 's', 'MarkerFaceColor', h2.CData); hold on;

% PDW phase
Usd = [3, 4, 6, 8, 2, 4, 6, 8, 8];
delta = [1/32, 1/32, 1/32, 1/32, 1/8, 1/8, 1/8, 1/8, 1/2];
h3 = scatter(Usd, delta,my_marker_size, "filled"); hold on;

Usd = [8, 2];
delta = [ 1-0.6136822006913542, 1-0.7038300868950937];
scatter(Usd, delta,my_marker_size, "filled", 's', 'MarkerFaceColor', h3.CData); hold on;

% uniform-SC phase
Usd = [0, 0, 0.3, 1];
delta = [1/32, 1/8, 1/8, 1/8];
h4 = scatter(Usd, delta,my_marker_size, "filled");

Usd = [0.3];
delta = [1-0.5126124251571303];
scatter(Usd, delta, my_marker_size,"filled", 's', 'MarkerFaceColor', h4.CData); hold on;

set(gca, 'fontsize', 20);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2);
xlabel('$U_{sd}$', 'Interpreter', 'latex');
ylabel('$\delta$', 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'FontSize', 20);
set(get(gca, 'YLabel'), 'FontSize', 20);
box on;