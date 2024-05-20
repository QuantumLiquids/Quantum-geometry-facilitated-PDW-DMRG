U = 8;
colororder("gem");
C = colororder;

% Define marker size
my_marker_size = 100;

% fully gapped
Usd = [2,4,6,8];
delta = [0,0,0,0];
h1 = scatter(Usd, delta, my_marker_size, "filled"); hold on;

Usd = [8, 2];
nf = [1, 0.9999700312373689];
delta = 1 -nf;
scatter(Usd, delta, my_marker_size, "filled", 's', 'MarkerFaceColor', h1.CData); hold on;

% (pi,pi) CDW
Usd = [0, 1];
delta = [0,1/32];
h2 = scatter(Usd, delta, my_marker_size, "filled"); hold on;

Usd = [0.3, 0.3,1];
nf = [1, 0.9996862959558018, 0.9999038227242601];
delta = 1-nf;
% Use the same color as h2 for h2s and set marker size to my_marker_size
h2s = scatter(Usd, delta, my_marker_size, "filled", 's', 'MarkerFaceColor', h2.CData); hold on;

% PDW phase
Usd = [3, 4, 6, 8, 2, 4, 6, 8, 2,4,6 8];
delta = [1/32, 1/32, 1/32, 1/32, 1/8, 1/8, 1/8, 1/8, 1/2,1/2,1/2,1/2];
h3 = scatter(Usd, delta,my_marker_size, "filled"); hold on;

Usd = [8, 8,4,  2,2];
nf = [0.6141172698111564, 0.7857378284644683,  ...
    0.6477169003514064, ...
    0.8832484910611048,0.7038300868950937];
delta = 1-nf;
scatter(Usd, delta,my_marker_size, "filled", 's', 'MarkerFaceColor', h3.CData); hold on;

% uniform-SC phase
Usd = [0, 0, 0.3, 1];
delta = [1/32, 1/8, 1/8, 1/8];
h4 = scatter(Usd, delta,my_marker_size, "filled");

Usd = [0.3,0.3];
nf = [0.5492084164030157, 0.7218687596503331];
delta = 1-nf ;
scatter(Usd, delta, my_marker_size,"filled", 's', 'MarkerFaceColor', h4.CData); hold on;

set(gca, 'fontsize', 20);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2);
xlabel('$U_{sd}$', 'Interpreter', 'latex');
ylabel('$\delta$', 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'FontSize', 20);
set(get(gca, 'YLabel'), 'FontSize', 20);
box on;