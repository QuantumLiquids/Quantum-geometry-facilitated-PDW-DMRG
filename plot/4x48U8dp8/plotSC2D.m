Ly = 4;
Lx = 48;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 8;
Hole = Ly*Lx*2/8;
D_values = [7000,10000,12000,15000,18000];

% trunc_errs = [ 1.75e-07, 1.29e-07, 9.92e-08,6.5e-08,4.64e-08]';
trunc_errs = 1./D_values;
bond_width = 2;
bond_color = 'k';
positive_sc_color = [0     172	156]/256;
negative_sc_color = [237   159	84 ]/256;
reference_point_color = [152  115	232]/256;
circle_scale = 2;
legend_entries = cell(size(D_values));

% x_end = Lx*3/4+1;
x_end = Lx/4+14;
for i = 1:numel(D_values)
    if i == numel(D_values)
        D = D_values(i);
        % Create the file path
        file_path = ['../../data/onsitepair', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
            'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
            'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

        % Load the data from the JSON file
        data = jsondecode(fileread(file_path));

        % Filter the data based on data{i}{1}(1) == Lx * Ly / 2, the
        % reference site is the s-orbital
        sc_s = {};
        sc_d = {};
        counts = 1;
        countd = 1;
        for j = 1:numel(data)
            if data{j}{1}(1) == Lx * Ly / 2
                if mod(data{j}{1}(2) - data{j}{1}(1), 2 ) == 0
                    sc_s{counts} = data{j};
                    counts = counts + 1;
                else
                    sc_d{countd} = data{j};
                    countd = countd + 1;
                end
            end
        end
        % ===== s orbital ===== %
        figure;
        for x = Lx/4+1:x_end
            for y= 1:Ly-1
                line([x,x],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;  
                if x == Lx*3/4+1
                    line([x+1,x+1],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
                end
            end
            for y = 1:Ly
                line([x,x+1],[y,y],'color',bond_color,'linewidth',bond_width);  hold on;
            end
        end


        for j = 1:numel(sc_s)
            site_idx = sc_s{j}{1}(2); % C++ convention
            s_orbital_sc_correlation = sc_s{j}{2};
            x = fix(site_idx /(2*Ly)) + 1;
            if(x > x_end+1)
                continue
            end
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
        ref_site_idx = sc_s{1}{1}(1); % C++ convention
        x = fix(ref_site_idx /(2*Ly)) + 1;
        y = mod(ref_site_idx, 2 * Ly) /2 +1;
        radius = sqrt(abs(s_orbital_sc_correlation)) * 1.5 * circle_scale;
        side_length = 2* radius;
        rectangle('Position', [x - side_length/2, y - side_length/2, side_length, side_length],...
          'Curvature', 0, 'FaceColor', 'none', 'EdgeColor', reference_point_color, 'LineWidth', 3);

        axis off;
        axis equal;

        %====== d orbital ======%
        figure;
        for x = Lx/4+1:x_end
            for y= 1:Ly-1
                line([x,x],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;  
                if x == Lx*3/4+1
                    line([x+1,x+1],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
                end
            end
            for y = 1:Ly
                line([x,x+1],[y,y],'color',bond_color,'linewidth',bond_width);  hold on;
            end
        end

        for j = 1:numel(sc_d)
            site_idx = sc_d{j}{1}(2); % C++ convention
            d_orbital_sc_correlation = sc_d{j}{2};
            x = fix((site_idx -1) /(2*Ly)) + 1;
            if(x > x_end+1)
                continue
            end
            y = mod((site_idx -1), 2 * Ly) /2 +1;
            center = [x, y];
            radius = sqrt(abs(d_orbital_sc_correlation)) * circle_scale;
            if d_orbital_sc_correlation >= 0
                % Draw the disk using the rectangle function
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', positive_sc_color, 'EdgeColor', 'none');
            else
                rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius],...
                    'Curvature', [1, 1], 'FaceColor', negative_sc_color, 'EdgeColor', 'none');
            end
        end
        ref_site_idx = sc_d{1}{1}(1); % C++ convention
        x = fix(ref_site_idx /(2*Ly)) + 1;
        y = mod(ref_site_idx, 2 * Ly) /2 +1;
        radius = sqrt(abs(d_orbital_sc_correlation)) * 1.5 * circle_scale;
        side_length = 2 * radius;
        rectangle('Position', [x - side_length/2, y - side_length/2, side_length, side_length],...
          'Curvature', 0, 'FaceColor', 'none', 'EdgeColor', reference_point_color, 'LineWidth', 3);

        axis off;
        axis equal;
    end
end


