%% Run this code to obtain the result in Fig. 11.
clc
clear 
close all 
addpath(genpath('.\Comparison_methods\'));
addpath(genpath('.\Proposed_methods\'));
addpath(genpath('.\Utils\'));

N1 = 19; % size of the 2D spectrally sparse signal 
N2 = 19; % size of the 2D spectrally sparse signal 
r = 4; % rank
delta_f = 30*1e3; % subcarrier spacing
T = 1 / delta_f; % effective OFDM symbol duration
Tcp = 0.5 / delta_f; % cyclic prefix duration
Ts = T + Tcp; % total OFDM symbol duration
c0 = 3e+8;  % speed of light
fc = 20*1e9; % carrier frequency
lambda = c0 / fc; % wavelength
SNR_vec = [5,10,15,20,25]; % SNR vector

max_iter = 1000; % maximum number of iterations
tol = 1e-10; % tolerance 
m = 50; % number of samples      
p = m/(N1*N2); % sampling rate
num_trial = 100; % Monte Carlo trials

NRE1 = zeros(length(SNR_vec),num_trial); 
NRE2 = NRE1;
NRE3 = NRE1;
for i = 1:length(SNR_vec)
    SNR = SNR_vec(i);
    fprintf('SNR = %d dB\n', SNR);
    parfor trial = 1:num_trial
        rng(trial);
        Range = 200 * rand(r, 1); % range of the target         
        Veloc = 30 * rand(r, 1); % velocity of the target           
        Doppl = 2 * Veloc / lambda; % Doppler shift of the target
        Delay = 2 * Range.' / c0; % round-trip time-delay of the target  
        power = rand(r,1); % magnitude of target’s complex gain parameter
        phase = rand(r,1)*360; % phase of target’s complex gain parameter
        amp = 10.^(power/10).*exp(1i*phase*pi/180); % It corresponds to the parameter 'alpha_l' in the paper.
        
        [~,K,x_star] = generate_observation_2D(m,[N1 N2],r,amp,delta_f,Ts,Doppl,Delay); % generate the signal
        [x_star_n,~] = add_noise_to_echo(x_star,trial,SNR); % add the noise
        xs_n = zeros(size(x_star_n));
        xs_n(K) = x_star_n(K);
        X_star = reshape(x_star,[N1 N2]);
        Xs_n = reshape(xs_n,[N1 N2]);
        Obs_Hankel = H2_2D(Xs_n);
        omega = find(Obs_Hankel);
        Hx_GT = H2_2D(X_star);
        
        %% SHGD method
        opt = 0;
        stepsize = 0.5;
        [X1 ,~,~,~] = SHGD_2D(Xs_n(1:N1,1:N2),K,N1,N2,r,p,tol,max_iter,opt,stepsize,X_star);
        
        %% SUMC and SUMC-L methods
        rho = 1.05;
        mu = 1e-4;
        w1 = 0.6;
        w2 = 1;
        ITERS_num = 5;
        thr = 0.9;
        eps = norm(Hx_GT(omega)-Obs_Hankel(omega));
        knowledge = init_2D(Xs_n(1:N1,1:N2),r); % initialization
        X2 = SUMC_2D(Xs_n(1:N1,1:N2),knowledge,max_iter,rho,mu,tol,r,w1,ITERS_num,thr*eps);
        X3 = SUMC_L_2D(Xs_n(1:N1,1:N2),knowledge,max_iter,rho,mu,tol,r,w2,ITERS_num,thr*eps);
        
        NRE1(i,trial) = norm((X1(1:N1,1:N2) - X_star),'fro')/norm(X_star,'fro');
        NRE2(i,trial) = norm((X2(1:N1,1:N2) - X_star),'fro')/norm(X_star,'fro');
        NRE3(i,trial) = norm((X3(1:N1,1:N2) - X_star),'fro')/norm(X_star,'fro');
    end
end

figure(11); % Fig. 11 in the paper
semilogy(SNR_vec,mean(NRE1,2),'LineWidth',3,'Marker','^','MarkerSize',9,'color',[0.93,0.69,0.13]); hold on
semilogy(SNR_vec,mean(NRE2,2),'b','LineWidth',3,'Marker','o','MarkerSize',9); hold on
semilogy(SNR_vec,mean(NRE3,2),'r','LineWidth',3,'Marker','square','MarkerSize',9); hold on
xlabel('SNR (dB)','fontsize',17,'FontName','Times new roman');
ylabel('NRE','fontsize',17,'FontName','Times new roman');
grid on
set(gca, 'FontName', 'Times new roman', 'FontSize', 17);
h = legend('SHGD','SUMC','SUMC-L');
set(h,'Interpreter','latex','fontsize',14.5,'FontName','Times new roman','Position',[0.63036 0.69286 0.2370 0.1991]);