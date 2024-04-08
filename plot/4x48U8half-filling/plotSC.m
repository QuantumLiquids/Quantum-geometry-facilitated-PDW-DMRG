Ly = 4;
Lx = 48;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 8;
Hole = 0;
D_values = [5000,7000,10000];

% trunc_errs = [ 4.27e-08,2.72e-08, 1.49e-08]';
sc_corr_finite_D = [];
fit_length = 20;
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
    semilogy(x_values, y_values, 'x', 'MarkerSize', 6);
    hold on;
    
    sc_corr_finite_D = [sc_corr_finite_D; y_values];
    % Generate the legend entry for the current D value
    if i == 1
        legend_entries{i} = ['$D = ', num2str(D),'$'];
    else
        legend_entries{i} = ['$', num2str(D),'$'];
    end
end

sc_extraplt = sc_corr_finite_D(end,:);

% Fit a exponential function to SC correlation
log_sc_extraplt = log(sc_extraplt);

% Perform linear regression
X = [ones(length(x_values), 1), x_values']; % Design matrix
coefficients = X \ log_sc_extraplt'; % Coefficients of the linear model

% Extract fitted parameters
intercept = coefficients(1);
slope = coefficients(2);
xi = -1/slope;
fprintf('Fitted parameters:\n');
fprintf('correlation length xi : %f\n', xi);
% Generate fitted curve
fitted_curve = X * coefficients;

plot(x_values, exp(fitted_curve), '-.');
hold off;

% Set the labels and title
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$r$','Interpreter','latex');
%Phi(x) = \langle\Delta(0)^\dagger \Delta(x)\rangle
ylabel('$\Phi(r)$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% Display the legend
% l=legend(legend_entries, 'Location', 'best');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');

% ylim([1e-3, 1e-1]);
xlim([0 25])
% xticks([2,4,8,16,24]);
% Display the plot
% grid on;