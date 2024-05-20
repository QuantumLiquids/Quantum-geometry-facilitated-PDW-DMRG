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
D=20000;

marker_size = 8;

marker_colors{1} = [019, 103, 131]/256;
marker_colors{2} = [255,158,002] / 256;
marker_colors{3} = [251,056,071] / 256;
marker_colors{4} = [131,064,028] / 256;
marker_colors{5} = [075,116,178] / 256;
marker_colors{6} = [107,112,092] / 256;

%SC
file_path = ['../../data/onsitepair', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
    'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
    'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];
data = jsondecode(fileread(file_path));
% Filter the data based on data{i}{1}(1) == Lx * Ly / 2
filtered_scss_data = {};
filtered_scsd_data = {};
count = 1;
for j = 1:numel(data)
    if data{j}{1}(1) == Lx * Ly / 2 && mod(data{j}{1}(2) - data{j}{1}(1), 2 * Ly) == 0
        filtered_scss_data{count} = data{j};
    end
    if data{j}{1}(1) == Lx * Ly / 2 && mod(data{j}{1}(2) - data{j}{1}(1), 2 * Ly) == 1
        filtered_scsd_data{count} = data{j};
        count = count + 1;
    end
end

% Extract the x and y values using a for loop
x_values = zeros(1, numel(filtered_scss_data));
sc_ss = zeros(1, numel(filtered_scss_data));
sc_sd = zeros(1, numel(filtered_scsd_data));
for j = 1:numel(filtered_scss_data)
    x_values(j) = (filtered_scss_data{j}{1}(2) - filtered_scss_data{j}{1}(1)) / (2*Ly);
    sc_ss(j) = filtered_scss_data{j}{2};
end

for j = 1:numel(filtered_scsd_data)
    sc_sd(j) = filtered_scsd_data{j}{2};
end


% Compute the logarithm of the data
log_sc_ss = log(sc_ss);
log_sc_sd = log(-sc_sd);

% Perform linear fits
p_sc_ss = polyfit(x_values, log_sc_ss, 1);
p_sc_sd = polyfit(x_values, log_sc_sd, 1);

% Extract the parameters
a_sc_ss = exp(p_sc_ss(2)); % Intercept
xi_sc = -1 / p_sc_ss(1); % Slope

a_sc_sd = exp(p_sc_sd(2)); % Intercept
lambda_sc_sd = -1 / p_sc_sd(1); % Slope

% Generate fitted data
fitted_ss = a_sc_ss * exp(-x_values / xi_sc);
fitted_sd = a_sc_sd * exp(-x_values / lambda_sc_sd);

% Plot the fitted data
h_fit_ss = semilogy(x_values, fitted_ss, '--', 'LineWidth', 1.5,'Color',marker_colors{6}); hold on;

% Plot the data on a logarithmic scale
h_sc_ss = semilogy(x_values, sc_ss, 'o', 'MarkerSize', marker_size);hold on;
h_sc_sd = semilogy(x_values, -sc_sd, 'x', 'MarkerSize', marker_size);


%Charge density correlation

% Create the file path
file_path = ['../../data/nfnf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
    'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
    'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

% Load the data from the JSON file
corr_data = jsondecode(fileread(file_path));

file_path = ['../../data/nf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
    'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
    'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

% Load the data from the JSON file
nf_data = jsondecode(fileread(file_path));

markers_band = ['s','+'];
band_name = {'s-orbital','d-orbital'};
for band = [0,1] %s and d orbital
    % Filter the data based on data{i}{1}(1) == Lx * Ly / 2
    filtered_data = {};
    count = 1;
    for j = 1:numel(corr_data)
        site1 = corr_data{j}{1}(1);
        site2 = corr_data{j}{1}(2);
        if site1 == Lx * Ly / 2 +band && mod(site2-site1, 2 * Ly) == 0
            filtered_data{count} = corr_data{j};
            count = count + 1;
        end
    end
    trunc_length = 5;
    % the distance and the connected correlation
    x_values = zeros(1, trunc_length);
    charge_correlation_values = zeros(1, trunc_length);
    for j = 1:trunc_length
        site1 = filtered_data{j}{1}(1);
        site2 = filtered_data{j}{1}(2);
        x_values(j) = (filtered_data{j}{1}(2) - filtered_data{j}{1}(1)) / (2*Ly);
        charge_correlation_values(j) = filtered_data{j}{2} - nf_data(site1 + 1, 2) * nf_data(site2 + 1,2);
    end

    % Plot the data on a logarithmic scale
    h(band+1) = semilogy(x_values, abs(charge_correlation_values), markers_band(band+1), 'MarkerSize', marker_size);
end
h_charge_ss = h(1);
h_charge_sd = h(2);


%fitting
% Compute the logarithm of the data
log_correlation = log(abs(charge_correlation_values));

% Perform linear fits
p_charge = polyfit(x_values, log_correlation, 1);

% Extract the parameters
A = exp(p_charge(2)); % Intercept
xi_charge = -1 / p_charge(1); % Slope

% Generate fitted data
fitted_data = A * exp(-[2:1:20]/ xi_charge);

% Plot the fitted data
h_fit_charge = semilogy([2:1:20], fitted_data, '--', 'LineWidth', 1.5,'Color',marker_colors{6}); hold on;
uistack(h_fit_charge,'bottom')

%spin

file_path1 = ['../../data/zzsf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
    'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
    'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];
% Load the data from the JSON file
data1 = jsondecode(fileread(file_path1));

file_path2 = ['../../data/pmsf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
    'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
    'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];
% Load the data from the JSON file
data2 = jsondecode(fileread(file_path2));

file_path3 = ['../../data/mpsf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
    'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
    'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];
% Load the data from the JSON file
data3 = jsondecode(fileread(file_path3));

% Filter the data
filtered_datas = {};
filtered_datad = {};
count = 1;
for ii = 1:length(data1)
    entry = data1{ii};
    entry2 = data2{ii};
    entry3 = data3{ii};
    if (entry{1}(1) == Lx * Ly / 2 && mod(entry{1}(2) - entry{1}(1), 2 * Ly) == 0)
        entry{2} =  entry{2}+ 0.5 * (entry2{2} + entry3{2});
        filtered_datas{count} = entry;
    end

    if (entry{1}(1) == Lx * Ly / 2 && mod(entry{1}(2) - entry{1}(1), 2 * Ly) == 1)
        entry{2} =  entry{2}+ 0.5 * (entry2{2} + entry3{2});
        filtered_datad{count} = entry ;
        count = count + 1;
    end
end

% Extract the x and y values
x_value = (cellfun(@(entry) (entry{1}(2) - entry{1}(1)) / (2 * Ly), filtered_datas)).';
spin_correlation_ss = -cellfun(@(entry) entry{2}, filtered_datas).';
spin_correlation_sd = -cellfun(@(entry) entry{2}, filtered_datad).';


h_spin_ss = semilogy(x_value, abs(spin_correlation_ss),'<', 'MarkerSize',marker_size);
h_spin_sd = semilogy(x_value, abs(spin_correlation_sd),'>','MarkerSize',marker_size);
%fitting
% Compute the logarithm of the data
log_correlation = log(abs(spin_correlation_ss(1:10)));
p_spin_ss = polyfit(x_value(1:10), log_correlation, 1);
A = exp(p_spin_ss(2)); % Intercept
xi_s1 = -1 / p_spin_ss(1); % Slope
fitted_data = A * exp(-[2:1:20]/ xi_s1);
h_fit_spin = semilogy([2:1:20], fitted_data, '--', 'LineWidth', 1.5,'Color',marker_colors{6}); hold on;
uistack(h_fit_spin,'bottom')

log_correlation = log(abs(spin_correlation_sd(1:10)));
p_spin_ss = polyfit(x_value(1:10), log_correlation, 1);
A = exp(p_spin_ss(2)); % Intercept
xi_s2 = -1 / p_spin_ss(1); % Slope
fitted_data = A * exp(-[2:1:20]/ xi_s2);
h_fit_spin = semilogy([2:1:20], fitted_data, '--', 'LineWidth', 1.5,'Color',marker_colors{6}); hold on;
uistack(h_fit_spin,'bottom')


ylim([1e-15,1e-0]);
xlim([2,12]);

set(gca,'fontsize',20);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$ x$','Interpreter','latex');
ylabel('Correlations','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',20);
set(get(gca,'YLabel'),'FontSize',20);

l = legend([h_sc_ss,  h_charge_ss, h_spin_ss,h_sc_sd, h_charge_sd,h_spin_sd], ...
    {'$\Phi^{ss}$', '$|D^{ss}|$', '$|F^{ss}|$',  ...
    ['$-\Phi^{sd}\quad \xi_{SC} = ',num2str(xi_sc,'%.2f'),'$'], ...
    ['$|D^{sd}|\quad \xi_{C} = ',num2str(xi_charge,'%.2f'),'$'], ...
    ['$|F^{sd}|\quad \xi_{S} = ',num2str(xi_s1,'%.2f'),',' ,num2str(xi_s2,'%.2f'),'$']}, ...
   'NumColumns', 2);

set(l,'Location', 'best');
set(l, 'Box', 'off');
set(l, 'Interpreter', 'latex');
set(l, 'FontSize', 18);


