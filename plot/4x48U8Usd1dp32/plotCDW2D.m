Ly = 4;
Lx = 48;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 1;
Hole = Lx * Ly * 2/32;
D_values = [10000];


bond_width = 2;
bond_color = 'k';
% color1 = [233, 196, 107]/256;
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

        figure;
        for x = 1:Lx
            for y= 1:Ly
                if y < Ly
                    line([x,x],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
                end
                if x < Lx
                    line([x,x+1],[y,y],'color',bond_color,'linewidth',bond_width);  hold on;
                end
                center = [x, y];
                radius = abs(s_density(y,x)) * circle_scale;
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', color1, 'EdgeColor', 'none');
            end
        end

        axis off;
        axis equal;

        %====== d orbital ======%
        figure;

        for x = 1:Lx
            for y= 1:Ly
                if y < Ly
                    line([x,x],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
                end
                if x < Lx
                    line([x,x+1],[y,y],'color',bond_color,'linewidth',bond_width);  hold on;
                end
                center = [x, y];
                radius = abs(d_density(y,x)) * circle_scale;
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', color1, 'EdgeColor', 'none');
            end
        end
        axis off;
        axis equal;
    end
end
