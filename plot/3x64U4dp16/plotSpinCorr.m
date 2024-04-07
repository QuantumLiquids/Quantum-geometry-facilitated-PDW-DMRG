clear all;
figure;
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
D_values = [10000,12000,14000,16000,20000];
trunc_errs = [2.51e-07, 1.75e-07, 1.29e-07, 9.92e-08,6.5e-08]';

spin_corr_finite_D_s = [];
spin_corr_finite_D_d = [];
end_point = 16;
for idx = 1:length(D_values)
    D = D_values(idx);
    % Create the file path
    file_path1 = ['../../data/zzsf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd, '%.1f'), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];
    % Load the data from the JSON file
    data1 = jsondecode(fileread(file_path1));

    file_path2 = ['../../data/pmsf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd, '%.1f'), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];
    % Load the data from the JSON file
    data2 = jsondecode(fileread(file_path2));

    file_path3 = ['../../data/mpsf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd, '%.1f'), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];
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
    y_values = -cellfun(@(entry) entry{2}, filtered_datas).';
    y_valued = -cellfun(@(entry) entry{2}, filtered_datad).';
    
    x_value = x_value(1:end_point);
    y_values = y_values(1:end_point)';
    y_valued = y_valued(1:end_point)';
    spin_corr_finite_D_s = [spin_corr_finite_D_s; y_values];
    spin_corr_finite_D_d = [spin_corr_finite_D_d; y_valued];
    % Plot the data on a logarithmic scale
    % semilogy(x_values, y_values, 'o','MarkerSize', 8);
    % hold on;
end

spin_extraplts = zeros(1, size(spin_corr_finite_D_s, 2));
spin_extrapltd = zeros(1, size(spin_corr_finite_D_d, 2));

for col = 1:size(spin_corr_finite_D_s, 2)
    p = polyfit(trunc_errs, spin_corr_finite_D_s(:, col), 2);
    spin_extraplts(col) = polyval(p, 0);

    p = polyfit(trunc_errs, spin_corr_finite_D_d(:, col), 2);
    spin_extrapltd(col) = polyval(p, 0);
end
semilogy(x_value, spin_extraplts, 'o', 'MarkerSize', 8); hold on;
% semilogy(x_value, abs(spin_extrapltd), 'x', 'MarkerSize', 8); hold on;

% exponentially decay fit
% Take the logarithm of the spin correlation function
log_spin_extraplt = log(spin_extraplts);

% Perform linear regression
X = [ones(length(x_value), 1), x_value]; % Design matrix
coefficients = X \ log_spin_extraplt'; % Coefficients of the linear model

% Extract fitted parameters
intercept = coefficients(1);
slope = coefficients(2);
xi = -1/slope;
fprintf('correlation length xi: %.4f\n', xi);
% Generate fitted curve
fitted_curve = X * coefficients;

plot(x_value, exp(fitted_curve), '-.');

ylim([1e-15, 1e-3]);
xlim([2 32]);

hold off;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$r$','Interpreter','latex');

ylabel('$-F(r)$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
