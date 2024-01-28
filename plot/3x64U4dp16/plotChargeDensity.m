% Define the initial values
Ly = 3;
Lx = 64;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 3.8;
Udd = 3.7;
Usd = 4.0;
Hole = 24;
D_values = [6000];

for i = 1:length(D_values)
    D = D_values(i);
    % Create the file path
    file_path = ['../data/nf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd, '%.1f'), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    data = jsondecode(fileread(file_path));

    s_data = data(mod(data(:, 1), 2*Ly) == 0, :);
    d_data = data(mod(data(:, 1), 2*Ly) == 1, :);
    
    % Extract the x and y values
    x_values = s_data(:, 1) / (2*Ly);
    n_s = s_data(:, 2);
    n_d = d_data(:, 2);
    
    plot(x_values, n_s, '-o', 'DisplayName', ['$n_s, D = ', num2str(D),'$']);
    hold on;
    plot(x_values, n_d, '-x', 'DisplayName', ['$n_d, D = ', num2str(D),'$']);
end

hold off;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$ x$','Interpreter','latex');
ylabel('$\langle n(x)\rangle$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);


l=legend('Location', 'best');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');
%Display the plot
%grid on;