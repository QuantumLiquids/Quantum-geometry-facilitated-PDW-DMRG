Ly = 4;
Lx = 36;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 0;
Hole = 0;
D_values = [20000];

bond_width = 2;
bond_color = 'k';
color1 = [042, 157, 142]/256;
circle_scale = 0.3;
legend_entries = cell(size(D_values));

for i = 1:numel(D_values)
    if i == numel(D_values)
        D = D_values(i);
        % Create the file path
        file_path = ['../../data/nf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
            'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
            'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

        % Load the data from the JSON file
        data = jsondecode(fileread(file_path));

        s_data = data(mod(data(:, 1), 2) == 0, :);
        s_density = reshape(s_data(:,2), Ly,[]);
        d_data = data(mod(data(:, 1), 2) == 1, :);
        d_density = reshape(d_data(:,2), Ly,[]);

        % Determine the global color limits
        combined_data = [s_density(:); d_density(:)];
        clim = [min(combined_data), max(combined_data)];

        % Create a new figure with subplots
        figure;

        % Subplot for s_density
        subplot(2, 1, 1);
        imagesc(1:Lx, 1:Ly, s_density, clim);
        title('s-orbital density');
        set(gca, 'YDir', 'normal');
        axis equal;
        axis off;

        % Subplot for d_density
        subplot(2, 1, 2);
        imagesc(1:Lx, 1:Ly, d_density, clim);
        title('d-orbital density');
        set(gca, 'YDir', 'normal');
        axis equal;
        axis off;

        % Create a single colorbar for both subplots
        % Position it to the right side of the figure
        h = colorbar;
        colormap(flipud(colormap('summer')));
        h.Ticks = [0.6, 0.8,1.0,1.2,1.4];
        h.FontName ='Arial';
        h.FontSize = 14;
        set(h, 'Position', [0.92 0.11 0.02 0.815]);
        

        % Adjust the position of the colorbar to align with the subplots
        pos1 = get(subplot(2, 1, 1), 'Position');
        pos2 = get(subplot(2, 1, 2), 'Position');
        set(subplot(2, 1, 1), 'Position', [pos1(1) pos1(2) pos1(3) pos1(4) * 0.9]);
        set(subplot(2, 1, 2), 'Position', [pos2(1) pos2(2) pos2(3) pos2(4) * 0.9]);
    end
end
