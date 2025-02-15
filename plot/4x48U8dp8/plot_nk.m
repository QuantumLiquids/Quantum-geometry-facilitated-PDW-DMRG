
ky_idx = 2;
kz_idx = 2;
% Load nk data
% nkD5000 = jsondecode(fileread('nkD5000.json'));
nkD15000 = jsondecode(fileread('nkD15000.json'));


% kx_set for x-axis
kx_set = [-pi:0.1:pi,0.0];
kx_set = sort(kx_set);

% Plot the two curves
figure;

% plot(kx_set, nkD15000(:,ky_idx,kz_idx), '-', 'DisplayName', 'D = 15000');
plot(kx_set, nkD15000(:,1,1), '-', 'DisplayName', 'k_y=0, k_z=0'); hold on;
plot(kx_set, nkD15000(:,1,2), '-.', 'DisplayName', 'k_y=0, k_z=\pi');
plot(kx_set, nkD15000(:,2,1), '--', 'DisplayName', 'k_y=\pi/2, k_z=0');
plot(kx_set, nkD15000(:,2,2), '-.', 'DisplayName', 'k_y=\pi/2, k_z=\pi');
plot(kx_set, nkD15000(:,3,1), '--', 'DisplayName', 'k_y=\pi, k_z=0');
plot(kx_set, nkD15000(:,3,2), '-.', 'DisplayName', 'k_y=\pi, k_z=\pi');
plot(kx_set, nkD15000(:,4,1), '--', 'DisplayName', 'k_y=3\pi/2, k_z=0');
plot(kx_set, nkD15000(:,4,2), '-.', 'DisplayName', 'k_y=3\pi/2, k_z=\pi');

set(gca,'fontsize',30);
set(gca,'linewidth',2.5);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
xlabel('$k_x$','Interpreter','latex');
ylabel('$n_k$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',30); 
set(get(gca,'YLabel'),'FontSize',30); 
set(gca,'Xlim', [0,pi]);
set(gca,'Ylim', [0,2]);
set(gcf,'position',[1000,1000,500,500]);

legend show;
