% Data
xi = [6.069436, 6.642056, 7.440832, 8.109383, 9.135854];
D = [7000, 9000, 12000, 15000, 20000];

% Plot the data
plot(D, xi, 'o','MarkerSize',8);
hold on;

% Fit the data with a power-law relationship xi = a * D^alpha
p = polyfit(log(D), log(xi), 1);
alpha = p(1);
a = exp(p(2));

% Generate the fitted line
D_fit = linspace(min(D), max(D), 100);
xi_fit = a * D_fit.^alpha;
disp(['Alpha value: ', num2str(alpha)]);


% Plot the fitted line
plot(D_fit, xi_fit, '--', 'Color', [0.5, 0.5, 0.5]);

% Set the labels and title
set(gca, 'fontsize', 24);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2);
xlabel('bond dimension D');
ylabel('correlation length \xi');
set(get(gca, 'XLabel'), 'FontSize', 24);
set(get(gca, 'YLabel'), 'FontSize', 24);
xlim([0, inf]);
ylim([0, inf]);
