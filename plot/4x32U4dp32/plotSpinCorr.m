% Define the initial values
Ly = 3;
Lx = 24;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 3.8;
Udd = 3.7;
Usd = 4.0;
Hole = 8;
D_values = [10000, 12000, 14000, 16000,20000];

for idx = 1:length(D_values)
    D = D_values(idx);
    
    % Create the file path
     file_path = ['../../data/zzsf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd, '%.1f'), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    data = jsondecode(fileread(file_path));
    
    % Filter the data
    filtered_data = {};
    count = 1;
    for ii = 1:length(data)
        entry = data{ii};
        if (entry{1}(1) == Lx * Ly / 2 && mod(entry{1}(2) - entry{1}(1), 2 * Ly) == 0)
            filtered_data{count} = entry;
                count = count + 1;
        end
    end
    
    % Extract the x and y values
    x_values = (cellfun(@(entry) (entry{1}(2) - entry{1}(1)) / (2 * Ly), filtered_data)).';
    y_values = -cellfun(@(entry) entry{2}, filtered_data).';
    
    % Plot the data on a logarithmic scale
    semilogy(x_values, y_values, 'o');
    hold on;
end

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$-\langle S^z(0) S^z(x) \rangle$', 'Interpreter', 'latex');
% ylim([1e-5, 1e-1]);
hold off;
