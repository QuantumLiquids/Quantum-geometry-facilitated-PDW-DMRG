% Define the initial values
Ly = 3;
Lx = 64;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 8;
Hole = Lx * Ly * 2/8;
D_values = [5000,7000,9000,12000,15000];
% trunc_errs = [2.51e-07, 1.75e-07, 1.29e-07, 9.92e-08,6.5e-08,4.64e-08];
ns_finite_D_data = [];
nd_finite_D_data = [];
figure; % Create a new figure

for i = 1:length(D_values)
    D = D_values(i);
    % Create the file path
    file_path = ['../../data/nf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    data = jsondecode(fileread(file_path));

    s_data = data(mod(data(:, 1), 2*Ly) == 0, :);
    d_data = data(mod(data(:, 1), 2*Ly) == 1, :);
    
    % Extract the x and y values
    x_coor = s_data(:, 1) / (2*Ly);
    n_s = s_data(:, 2);
    n_s = (n_s + flip(n_s))/2;
    ns_finite_D_data = [ns_finite_D_data; n_s'];
    n_d = d_data(:, 2);
    n_d = (n_d + flip(n_d))/2;
    nd_finite_D_data = [nd_finite_D_data; n_d'];
    
    plot(x_coor+1, n_s, '-', 'DisplayName', ['$n_s, D = ', num2str(D),'$']); hold on;
    % plot(x_coor, n_d, 'x', 'DisplayName', ['$n_d, D = ', num2str(D),'$']);
end

hold off;

set(gca,'fontsize',28);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$ x$','Interpreter','latex');
ylabel('$\langle n(x)\rangle$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% xlim([0 Lx/2]);
ave_density_line = yline(1-1/9,'--');
ave_density_line.LineWidth = 2;

% Create the legend
l = legend( {'$D=5000$', '$7000$','$9000$', '$12000$', '$15000$'}, ...
            'Location', 'best', 'NumColumns', 2);

set(l,'Box','off');
set(l,'Interpreter','latex');
set(l,'FontSize',24);

set(gcf,'position',[1000,1000,450,350]);

% Display the plot
% grid on;