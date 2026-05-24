%% Run this code to obtain the results in Figs. 3-5.
clc
clear
close all
addpath(genpath('.\Comparison_methods\'));
addpath(genpath('.\Proposed_methods\'));
addpath(genpath('.\Utils\'));

n1 = 50; % size of the Hankel matrix
n2 = n1; % size of the Hankel matrix
n = n1 + n2 - 1; % length of the spectrally sparse signal 
vec_m = round(linspace(n*0.1,n*0.7,16)); % number of samples
vec_rate = vec_m/n; % sampling rate
vec_r = 2:2:30; % rank
Monte = 100; % Monte Carlo trials

idx = hankel(1:n1,n1:(n1-1)+n1);
Idx = idx(:);

num_method = 8;
RMSE_ang = zeros(length(vec_m),length(vec_r),Monte,num_method);
RMSE_amp = RMSE_ang;
NRE = RMSE_ang;

for i = 1:length(vec_m)
  m = vec_m(i);
  rate = vec_rate(i);
  for j = 1:length(vec_r)
    fprintf('(%d,%d)\n',i,j);
    r = vec_r(j);  
    parfor k = 1:Monte
        [xs,omega,x_star,freq,amp] = generate_observation_1D(n,m,r,k);
        x_est = zeros(n,num_method);
        
        %% OMP method
        grid_num = n*8;
        [A,hatA] = prepare_OMP(omega,grid_num,n);
        x_est(:,1) = OMP(xs,A,hatA,r);
        
        %% EMaC method
        rho = 1.05;
        mu = 1e-4;
        tol = 1e-10;
        max_iter = 1000;
        x_est(:,2) = EMAC(xs,n,omega,max_iter,rho,mu,tol,Idx);   
        
        %% PGD and SHGD methods
        opt = 1;
        stepsize = 0.5;
        x_est(:,3) = PGD(xs,omega,n1,n2,n,r,rate,tol,max_iter,opt,stepsize,x_star);
        x_est(:,4) = SHGD(xs,omega,n,r,rate,tol,max_iter,opt,stepsize,x_star);
        
        %% RSMC, SUMC, RSMC-L, and SUMC-L methods
        weight = 1;
        iter = 5;
        knowledge = init_1D(xs,n,r,Idx); % initialization
        x_est(:,5) = RSMC(xs,knowledge,n,omega,max_iter,rho,mu,tol,r,weight,Idx);  
        x_est(:,6) = SUMC(xs,knowledge,max_iter,n,omega,rho,mu,tol,r,weight,iter,Idx); 
        x_est(:,7) = RSMC_L(xs,knowledge,n,omega,max_iter,1,tol,r,weight,Idx); 
        x_est(:,8) = SUMC_L(xs,knowledge,max_iter,n,omega,rho,mu,tol,r,weight,iter,Idx);

        for p = 1:num_method
            [RMSE_ang(i,j,k,p),RMSE_amp(i,j,k,p),NRE(i,j,k,p)] = evaluation(x_est(:,p),x_star,freq,amp,n,r);
        end
    end
  end
end

close all
titles = [{'(a) OMP'},'(b) EMaC','(c) PGD','(d) SHGD','(e) RSMC','(f) SUMC','(g) RSMC-L','(h) SUMC-L'];
thr = 1e-4;
imgshows(3,0,1,NRE<thr,titles); % Fig. 3 in the paper
imgshows(4,0,1,RMSE_ang<thr,titles); % Fig. 4 in the paper
imgshows(5,0,1,RMSE_amp<thr,titles); % Fig. 5 in the paper