% Define the initial values
Ly = 3;
Lx = 64;
ts = 1;
td = -1;
tsd_xy = 1;
tsd_nn = 0;
Uss = 8;
Udd = 8;
Usd = 8;
Hole = Lx * Ly * 2/8;
D_values = [5000];
% trunc_errs = [2.51e-07, 1.75e-07, 1.29e-07, 9.92e-08,6.5e-08,4.64e-08];
ns_finite_D_data = [];
nd_finite_D_data = [];
for i = 1:length(D_values)
    D = D_values(i);
    % Create the file path
    file_path = ['../../data/nf', num2str(Ly), 'x', num2str(Lx), 'ts', num2str(ts), 'td', num2str(td), ...
        'tsd_xy', num2str(tsd_xy), 'tsd_nn', num2str(tsd_nn), 'Uss', num2str(Uss), 'Udd', num2str(Udd), ...
        'Usd', num2str(Usd), 'Hole', num2str(Hole), 'D', num2str(D), '.json'];

    % Load the data from the JSON file
    data = jsondecode(fileread(file_path));

    s_data = data(mod(data(:, 1), 2*Ly) == 0, :);
    d_data = data(mod(data(:, 1), 2*Ly) == 1, :);
    
    % Extract the x and y values
    x_coor = s_data(:, 1) / (2*Ly);
    n_s = s_data(:, 2);
    n_s = (n_s + flip(n_s))/2;
    ns_finite_D_data = [ns_finite_D_data; n_s'];
    n_d = d_data(:, 2);
    n_d = (n_d + flip(n_d))/2;
    nd_finite_D_data = [nd_finite_D_data; n_d'];
    
    plot(x_coor, n_s, 'o', 'DisplayName', ['$n_s, D = ', num2str(D),'$']); hold on;
    plot(x_coor, n_d, 'x', 'DisplayName', ['$n_d, D = ', num2str(D),'$']);
end


% Extrapolation
% ns_extraplt = zeros(1, size(ns_finite_D_data, 2));
% nd_extraplt = zeros(1, size(nd_finite_D_data, 2));
% for col = 1:size(ns_finite_D_data, 2)
%     p = polyfit(trunc_errs, ns_finite_D_data(:, col), 2);
%     ns_extraplt(col) = polyval(p, 0);
% 
%     p = polyfit(trunc_errs, nd_finite_D_data(:, col), 2);
%     nd_extraplt(col) = polyval(p, 0);
% end
% h1 = plot(x_coor+1, ns_extraplt, 'o', 'DisplayName', ['$n_s, D = ', num2str(D),'$'],'MarkerSize', 8); hold on;
% h2 = plot(x_coor+1, nd_extraplt, 'x', 'DisplayName', ['$n_d, D = ', num2str(D),'$'],'MarkerSize', 8);


% Fit Kc, plot Friedel oscillations
% n(x) - n_0 ~ A cos(Qx +phi) x^(-Kc/2)

% start_site = 3;
% fit_x = x_coor(start_site:Lx/2)+1;
% fit_ys = ns_extraplt(start_site : Lx/2);
% fit_yd = nd_extraplt(start_site : Lx/2);
% Q = 3*pi/8; %Momentum
% n0 = (Lx * Ly * 2 - Hole) / (Lx * Ly * 2);
% 
% model = @(params, x) params(1) * cos(params(5) * x + params(2)) .* x.^(-params(3)/2) + params(4);
% 
% % Initial guess for the parameters [A, phi, Kc]
% initial_guess = [1, 1, 1, n0, Q]; 
% 
% % s-orbital
% params_fit = lsqcurvefit(model, initial_guess, fit_x, fit_ys');
% A_cdw = params_fit(1);
% phi = params_fit(2);
% Kc = params_fit(3);
% x = fit_x(1):0.01:fit_x(end);
% fitted_curve = model(params_fit, x);
% 
% plot(x, fitted_curve, '-.'); hold on;
% % Display fitted parameters
% fprintf('Fitted parameters:\n');
% fprintf('Kc (s-orbital): %f\n', Kc);
% 
% % d-orbital
% params_fit = lsqcurvefit(model, initial_guess, fit_x, fit_yd');
% A_cdw = params_fit(1);
% phi = params_fit(2);
% Kc = params_fit(3);
% fitted_curve = model(params_fit, x);
% plot(x, fitted_curve, '-.'); hold on;
% % Display fitted parameters
% fprintf('Fitted parameters:\n');
% fprintf('Kc (d-orbital): %f\n', Kc);

hold off;

set(gca,'fontsize',28);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$ x$','Interpreter','latex');
ylabel('$\langle n(x)\rangle$','Interpreter','latex')
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% xlim([0 Lx/2]);

% l = legend([h1, h2], {'$s$-orbital', '$d$-orbital'}, 'Location', 'best');

% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');

%Display the plot
%grid on;