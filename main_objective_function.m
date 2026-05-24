%% Run this code to obtain the result in Fig. 1.
clc
clear
close all
addpath(genpath('.\Comparison_methods\'));
addpath(genpath('.\Proposed_methods\'));
addpath(genpath('.\Utils\'));

n1 = 2; % size of the Hankel matrix  
n2 = n1; % size of the Hankel matrix
n = n1+n2-1; % length of the spectrally sparse signal 
m = 2; % number of samples
r = 1; % rank 
[~,~,x_star,~,~] = generate_observation_1D(n,m,r,3); % x_star: ground-truth
omega = [1 3]'; %  index set of the observed entries
xs = zeros(size(x_star)); % incomplete observation
xs(omega) = x_star(omega);

real_min = -32;
real_max = 23;
imag_min = -23;
imag_max = 31;
grid1 = (abs(real_min)+real_max)*10+1;
grid2 = (abs(imag_min)+imag_max)*10+1;

Z = zeros(grid2,grid1);
real_part = linspace(real_min,real_max,grid1);
imag_part = linspace(imag_min,imag_max,grid2);
[X, Y] = meshgrid(real_part, imag_part);
w = 0.99; % weight hyperparameter
for i = 1:grid2
  for j = 1:grid1
    fprintf('(%d,%d)\n',i,j);
    x_per = x_star(2) + real_part(j)+1j*imag_part(i); % perturbation
    x_test = xs + x_per;
    x_test(omega) = xs(omega);
    [G, ~] = takagi(H2(x_test),r);
    M = (eye(n1)-w*G*G');
    Z(i,j) = norm_nuc(M*H2(x_test)*conj(M)); % objective function
  end
end

close all
figure(1); % Fig. 1 in the paper
set(gcf, 'Position', [100, 200, 550, 440]);
surf(X, Y, Z, 'edgecolor', 'none');
xlabel('${\rm{Re}}(x_2-y_2)$','interpreter','latex','fontsize',21,'FontName','Times new roman');
ylabel('${\rm{Im}}(x_2-y_2)$','interpreter','latex','fontsize',21,'FontName','Times new roman');
zlabel('$\|\mathcal{F}_w(\mathbf{x},\mathcal{E}_r(\mathbf{x}))\|_*$','interpreter','latex','fontsize',21,'FontName','Times new roman');
set(gca, 'FontSize', 21, 'FontName', 'Times new roman');
ax = gca;
ax.XLabel.Rotation = 9;
ax.YLabel.Rotation = -21;
ax.YLabel.Position = [-32.46,-0.11,-8.23];
xlim([real_min real_max]);
ylim([imag_min imag_max]);
set(gca, 'FontName', 'Times new roman', 'FontSize', 17, ...
    'XTick',[-30 0 20],'XTickLabel',{'-30','0','20'},...
    'YTick',[-20 0 30],'YTickLabel',{'-20','0','30'},...
    'ZTick',[0 20 40],'ZTickLabel',{'0','20','40'});
view([-35.5,18]);

