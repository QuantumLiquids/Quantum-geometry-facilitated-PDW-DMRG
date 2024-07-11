Ly = 4;
Lx = 24;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 0.3;
mu = -4.6;
D_values = [2000,4000,6000,8000,9000];

pairing_field_magnitudes = zeros(size(D_values));
for i = 1:length(D_values)
    D = D_values(i);
    % Create the file path
    file_path = ['../../data/Delta', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd), 'mus', num2str(mu),'mud', num2str(mu), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    data = jsondecode(fileread(file_path));
    pairing_field_magnitudes(i) = mean(abs(data(:,2)));
end

plot(D_values, pairing_field_magnitudes,'-x'); hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('bond dimension D');
ylabel('Pairing field')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

ylim([0,0.3]);