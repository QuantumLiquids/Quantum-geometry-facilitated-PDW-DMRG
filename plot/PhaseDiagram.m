U = 8;
colororder("gem");
C = colororder;

% Define marker size
my_marker_size = 60;
my_marker_size_square = 80;

% fully gapped
Usd = [2,4,6,8];
delta = [0,0,0,0];
h1 = scatter(Usd, delta, my_marker_size, "filled"); hold on;

Usd = [8, 2];
nf = [1, 0.9999700312373689];
delta = 1 -nf;
scatter(Usd, delta, my_marker_size_square, "filled", 's', 'MarkerFaceColor', h1.CData); hold on;

% (pi,pi) CDW
Usd = [0, 1];
delta = [0,1/32];
h2 = scatter(Usd, delta, my_marker_size, "filled"); hold on;

Usd = [0.3, 0.3,1,1];
nf = [1, 0.9996862959558018, 0.9999038227242601, 0.9330567123494481];
delta = 1-nf;
% Use the same color as h2 for h2s and set marker size to my_marker_size
h2s = scatter(Usd, delta, my_marker_size_square, "filled", 's', 'MarkerFaceColor', h2.CData); hold on;

% PDW phase
Usd = [3, 4, 6, 8, 2, 4, 6, 8, 6, 2,4,6 8];
delta = [1/32, 1/32, 1/32, 1/32, 1/8, 1/8, 1/8, 1/8, 1/4,1/2,1/2,1/2,1/2];
h3 = scatter(Usd, delta,my_marker_size, "filled"); hold on;

Usd = [8, 8,4,  2,2 , 1,1,1];
nf = [0.6141172698111564, 0.7857378284644683,  ...
    0.6477169003514064, ...
    0.8832484910611048,0.7038300868950937,...
    0.7598409960094165, 0.7029184804675991, 0.5253647279837241];
delta = 1-nf;
scatter(Usd, delta,my_marker_size_square, "filled", 's', 'MarkerFaceColor', h3.CData); hold on;

% uniform-SC phase
Usd = [0, 0, 0.3, 1,0.1];
delta = [1/32, 1/8, 1/8, 1/8,1/2];
h4 = scatter(Usd, delta,my_marker_size, "filled");

Usd = [0.3,0.3,0.3];
nf = [0.5648059912452188,  0.7455, 0.7819];
delta = 1-nf ;
scatter(Usd, delta, my_marker_size_square,"filled", 's', 'MarkerFaceColor', h4.CData); hold on;


%phase boundary line
%left - right line
x = [0.7, 0.6, 1, 1.3, 1.8];
y = [0.3, 0.5, 0.18, 0.125, 0];
y_fine = linspace(min(y), max(y), 100);

% Perform cubic spline interpolation
x_fine = spline(y, x, y_fine);

% Plot the smooth curve
plot(x_fine, y_fine, 'k-');
%PDW - gap line
x_start = spline(y, x, 0.015);
plot([x_start,8],[0.015, 0.015],'k-' );
%s-wave -CDW line

x_target = 1.3;
y_target = interp1(x_fine, y_fine, x_target, 'linear');
plot([0, x_target],[0.015, y_target],'k-' );

set(gca, 'fontsize', 20);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2);
xlabel('$U_{sd}$', 'Interpreter', 'latex');
ylabel('$\delta$', 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'FontSize', 20);
set(get(gca, 'YLabel'), 'FontSize', 20);
box on;