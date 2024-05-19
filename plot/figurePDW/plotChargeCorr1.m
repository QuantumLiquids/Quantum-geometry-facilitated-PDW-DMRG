
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
legend_entries = cell(size(D_values));
nnd_corr_finite_D = [];
nns_corr_finite_D = [];
for i = 1:numel(D_values)
    D = D_values(i);

    % Create the file path
    file_path = ['../../data/nfnf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd, '%.1f'), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    corr_data = jsondecode(fileread(file_path));

    file_path = ['../../data/nf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd, '%.1f'), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    nf_data = jsondecode(fileread(file_path));
    
    markers_band = ['o','x'];
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

        % the distance and the connected correlation
        x_values = zeros(1, numel(filtered_data));
        y_values = zeros(1, numel(filtered_data));
        for j = 1:numel(filtered_data)
            site1 = filtered_data{j}{1}(1);
            site2 = filtered_data{j}{1}(2);
            x_values(j) = (filtered_data{j}{1}(2) - filtered_data{j}{1}(1)) / (2*Ly);
            y_values(j) = filtered_data{j}{2} - nf_data(site1 + 1, 2) * nf_data(site2 + 1,2);
        end
        if(band == 0)
            nns_corr_finite_D = [nns_corr_finite_D;y_values];
        else
            nnd_corr_finite_D = [nnd_corr_finite_D;y_values];
        end
        % Plot the data on a logarithmic scale
        % loglog(x_values, abs(y_values), markers_band(band+1), 'MarkerSize', 6);
        % hold on;

        % Generate the legend entry for the current D value
        legend_entries{2 * i - 1 + band} = [band_name{band+1}, ', $D = ', num2str(D),'$' ];

    end
end


% Extrapolation
nns_extraplt = zeros(1, size(nns_corr_finite_D, 2));
nnd_extraplt = zeros(1, size(nnd_corr_finite_D, 2));

for col = 1:size(nnd_corr_finite_D, 2)
    p = polyfit(trunc_errs, nns_corr_finite_D(:, col), 2);
    nns_extraplt(col) = polyval(p, 0);
    p = polyfit(trunc_errs, nnd_corr_finite_D(:, col), 2);
    nnd_extraplt(col) = polyval(p, 0);
end
nns_extraplt = nns_extraplt * 10;
positiveIndices = nns_extraplt >= 0;
negativeIndices = nns_extraplt < 0;
h1 = loglog(x_values(positiveIndices), nns_extraplt(positiveIndices), 'o','MarkerSize', 6); hold on;
h1_negative = loglog(x_values(negativeIndices), -nns_extraplt(negativeIndices), 'o','MarkerSize', 6); hold on;

% h2 = loglog(x_values, abs(nnd_extraplt), 'x','MarkerSize', 8); hold on;

% Fit a power-law function to the last group of x_values and y_values
for band = 0
    if(band == 0)
        y_values =nns_extraplt;
    else
        y_values =nnd_extraplt;
    end
    indices = ismember(x_values, [3, 8, 9, 14,19,20,24,25]);
    log_x = log(x_values(indices));
    log_y = log(abs(y_values(indices)));
    fit = polyfit(log_x, log_y, 1);
    K = -fit(1);
    fprintf('Exponent Kc: %.4f\n', K);

    % Plot the fitted line
    x_guide = linspace(min(x_values), max(x_values), 100);
    y_guide = exp(polyval(fit, log(x_guide)));
    if(band == 0)
        h_fit1 = loglog(x_guide, y_guide, '--', 'LineWidth', 1.5);
    end
end




