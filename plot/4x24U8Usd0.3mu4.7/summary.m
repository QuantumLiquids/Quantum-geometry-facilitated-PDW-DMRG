color1 = [142 207 201]/256;
color2 = [255 190 122]/256;
color3 = [130 176 210]/256;
color4 = [019 033 060]/256;
color5 = [252 163 017]/256;
Ly = 4;
Lx = 24;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 0.3;
mu = -4.6;
D = 6000;

% charge density
file_path = ['../../data/nf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd), 'mus', num2str(mu),'mud', num2str(mu), 'D', num2str(D), '.json'];
% Load the data from the JSON file
data = jsondecode(fileread(file_path));

s_data = data(mod(data(:, 1), 2*Ly) == 0, :);
d_data = data(mod(data(:, 1), 2*Ly) == 1, :);

% Extract the x and y values
x_values = s_data(:, 1) / (2*Ly);
n_s = s_data(:, 2);
n_d = d_data(:, 2);
yyaxis left

plot(x_values+1, n_s, '-s', 'MarkerSize',8,'Color',color4);
hold on;

% SC
file_path = ['../../data/Delta', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
    'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
    'Usd', num2str(Usd), 'mus', num2str(mu),'mud', num2str(mu), 'D', num2str(D), '.json'];

ylabel('$n(x)$','Interpreter','latex')

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% Load the data from the JSON file
data = jsondecode(fileread(file_path));

s_data = data(mod(data(:, 1), 2*Ly) == 0, :);
d_data = data(mod(data(:, 1), 2*Ly) == 1, :);

% Extract the x and y values
x_values = s_data(:, 1) / (2*Ly);
Delta_s = s_data(:, 2);
Delta_d = d_data(:, 2);

yyaxis right
plot(x_values+1, Delta_s, '-o', 'MarkerSize',8,'Color',color5);
hold on;
plot(x_values+1, Delta_d, '-x', 'MarkerSize',8,'Color',color5);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$x$','Interpreter','latex');
ylabel('$\langle \Delta(x)\rangle$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% Match axis colors to line colors
ax = gca;
ax.YAxis(1).Color = color4; % Left y-axis color
ax.YAxis(2).Color = color5; % Right y-axis color

% l=legend('Location', 'best');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',16);
% set(l,'Location','SouthWest');