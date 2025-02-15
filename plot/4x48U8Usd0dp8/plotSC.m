Ly = 4;
Lx = 48;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 0;
Hole = Lx * Ly * 2/8;
D_values = [7000,9000,12000,15000,20000];

trunc_errs = 1./D_values;
sc_corr_finite_D = [];
fit_length = 25;
legend_entries = cell(size(D_values));

for i = 1:numel(D_values)
    D = D_values(i);

    % Create the file path
    file_path = ['../../data/onsitepair', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    data = jsondecode(fileread(file_path));

    % Filter the data based on data{i}{1}(1) == Lx * Ly / 2
    filtered_data = {};
    count = 1;
    for j = 1:numel(data)
        if data{j}{1}(1) == Lx * Ly / 2 && mod(data{j}{1}(2) - data{j}{1}(1), 2 * Ly) == 0
            filtered_data{count} = data{j};
            count = count + 1;
        end
    end

    % Extract the x and y values using a for loop
    x_values = zeros(1, numel(filtered_data));
    y_values = zeros(1, numel(filtered_data));
    for j = 1:numel(filtered_data)
        x_values(j) = (filtered_data{j}{1}(2) - filtered_data{j}{1}(1)) / (2*Ly);
        y_values(j) = filtered_data{j}{2};
    end

    % Plot the data on a logarithmic scale
    semilogy(x_values, y_values, '-o', 'MarkerSize', 6);
    hold on;

    % Generate the legend entry for the current D value
    legend_entries{i} = ['$D = ', num2str(D),'$'];

    % Fit a power-law function to the last group of x_values and y_values
    log_x = log(x_values(x_values<20));
    log_y = log(y_values(x_values<20));
    fit = polyfit(log_x, log_y, 1);
    K = -fit(1);
    fprintf('Exponent K: %.4f\n', K);

    % Plot the fitted line
    x_guide = linspace(min(x_values), max(x_values), 100);
    y_guide = exp(polyval(fit, log(x_guide)));

    % Fit a exponential function to SC correlation
    X = [ones(length(x_values(x_values<20)), 1), x_values(x_values<20)']; % Design matrix
    coefficients = X \ log_y'; % Coefficients of the linear model

    % Extract fitted parameters
    intercept = coefficients(1);
    slope = coefficients(2);
    xi = -1/slope;
    fprintf('correlation length xi : %f\n', xi);
   

    sc_corr_finite_D = [sc_corr_finite_D; y_values];
    % Generate the legend entry for the current D value
    if i == 1
        legend_entries{i} = ['$D = ', num2str(D),'$'];
    else
        legend_entries{i} = ['$', num2str(D),'$'];
    end
end


% Extrapolation
sc_extraplt = zeros(1, size(sc_corr_finite_D, 2));

for col = 1:size(sc_corr_finite_D, 2)
    p = polyfit(trunc_errs, sc_corr_finite_D(:, col), 2);
    sc_extraplt(col) = polyval(p, 0);
end
% loglog(x_values, sc_extraplt, '-o', 'MarkerSize', 8); hold on;

% Fit a power-law function to SC correlation

log_x = log(x_values(x_values<fit_length));
log_y = log(sc_extraplt(x_values<fit_length));
fit = polyfit(log_x, log_y, 1);
K = -fit(1);
fprintf('Exponent K: %.4f\n', K);

% Plot the fitted line
x_guide = linspace(min(x_values), max(x_values), 100);
y_guide = exp(polyval(fit, log(x_guide)));
% loglog(x_guide, y_guide, 'r--', 'LineWidth', 1.5);


% Set the labels and title
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$\Delta x$','Interpreter','latex');
ylabel('$\langle\Delta(0)^\dagger \Delta(x)\rangle$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% Display the legend
l=legend(legend_entries, 'Location', 'best');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');

% Display the plot
% grid on;