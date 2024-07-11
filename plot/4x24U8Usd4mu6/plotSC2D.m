Ly = 4;
Lx = 24;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 4;
mu = -6;
D_values = [10000];

bond_width = 2;
bond_color = 'k';
positive_sc_color = [0     172	156]/256;
negative_sc_color = [237   159	84 ]/256;
reference_point_color = [152  115	232]/256;
circle_scale = 1;
legend_entries = cell(size(D_values));

% x_end = Lx*3/4+1;
for i = 1:numel(D_values)
    if i == numel(D_values)
        D = D_values(i);
        % Create the file path
        file_path = ['../../data/Delta', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
            'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
            'Usd', num2str(Usd), 'mus', num2str(mu),'mud', num2str(mu), 'D', num2str(D), '.json'];

        % Load the data from the JSON file
        data = jsondecode(fileread(file_path));



        % Filter the data based on data{i}{1}(1) == Lx * Ly / 2, the
        % reference site is the s-orbital
        sc_s = data(1:2:end, :);
        sc_d = data(2:2:end, :);

        % ===== s orbital ===== %
        figure;
        for x = 1:Lx
            for y= 1:Ly-1
                line([x,x],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
            end
            for y = 1:Ly
                if x < Lx
                    line([x,x+1],[y,y],'color',bond_color,'linewidth',bond_width);  hold on;
                end
            end
        end


        for j = 1:size(sc_s, 1)
            site_idx = sc_s(j, 1); % C++ convention
            s_orbital_sc_correlation = sc_s(j,2);
            x = fix(site_idx /(2*Ly)) + 1;
            % if(x > x_end+1)
            %     continue
            % end
            y = mod(site_idx, 2 * Ly) /2 +1;
            center = [x, y];
            radius = sqrt(abs(s_orbital_sc_correlation)) * circle_scale;
            if s_orbital_sc_correlation >= 0
                % Draw the disk using the rectangle function
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', positive_sc_color, 'EdgeColor', 'none');
            else
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', negative_sc_color, 'EdgeColor', 'none');
            end
        end

        axis off;
        axis equal;

        %====== d orbital ======%
        figure;
        for x = 1:Lx
            for y= 1:Ly-1
                line([x,x],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
            end
            for y = 1:Ly
                if x < Lx
                    line([x,x+1],[y,y],'color',bond_color,'linewidth',bond_width);  hold on;
                end
            end
        end


        for j = 1:size(sc_d, 1)
            site_idx = sc_d(j, 1); % C++ convention
            s_orbital_sc_correlation = sc_d(j,2);
            x = fix(site_idx /(2*Ly)) + 1;
            % if(x > x_end+1)
            %     continue
            % end
            y = mod(site_idx, 2 * Ly) /2 +0.5;
            center = [x, y];
            radius = sqrt(abs(s_orbital_sc_correlation)) * circle_scale;
            if s_orbital_sc_correlation >= 0
                % Draw the disk using the rectangle function
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', positive_sc_color, 'EdgeColor', 'none');
            else
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', negative_sc_color, 'EdgeColor', 'none');
            end
        end

        axis off;
        axis equal;
    end
end


