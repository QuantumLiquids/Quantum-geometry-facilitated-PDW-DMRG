Ly = 4;
Lx = 48;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 8;
Hole = Lx * Ly * 2 / 32;
D_values = [7000,9000,12000,15000,18000,20000];
trunc_errs = [ 2.72e-08, 1.49e-08,1.50e-08,  9.68e-09, 7.27e-09,6.18e-09]' * 1e5;

corr_finite_D = [];
fit_length = 24;
legend_entries = cell(size(D_values));

for i = 1:numel(D_values)
    D = D_values(i);

    % Create the file path
    file_path = ['../../data/bupcbupa', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
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
    % semilogy(x_values, y_values, 'x', 'MarkerSize', 6); hold on;
    
    corr_finite_D = [corr_finite_D; y_values];
    % Generate the legend entry for the current D value
    if i == 1
        legend_entries{i} = ['$D = ', num2str(D),'$'];
    else
        legend_entries{i} = ['$', num2str(D),'$'];
    end
end


% Extrapolation
corr_extraplt = zeros(1, size(corr_finite_D, 2));

for col = 1:size(corr_finite_D, 2)
    p = polyfit(trunc_errs, corr_finite_D(:, col), 2);
    corr_extraplt(col) = polyval(p, 0);
end
h2 = semilogy(x_values, abs(corr_extraplt), '>', 'MarkerSize', 6); hold on;

% Fit a power-law function to correlation

x_value = (x_values(x_values<fit_length & mod(x_values, 2) ==0 ));
log_y = log(abs(corr_extraplt(x_values<fit_length & mod(x_values, 2) ==0)));

% Perform linear regression
X = [ones(length(x_value), 1), x_value']; % Design matrix
coefficients = X \ log_y'; % Coefficients of the linear model

% Extract fitted parameters
intercept = coefficients(1);
slope = coefficients(2);
xi = -1/slope;
fprintf('correlation length xi: %.4f\n', xi);
% Generate fitted curve
fitted_curve = X * coefficients;

h_fit2 = semilogy(x_value, exp(fitted_curve), '-.');
