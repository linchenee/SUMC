function [y,sigma_2] = add_noise_to_echo(x,trial,SNR)
% -------------------------------------------------------    
% This function is used to add the noise z to the echo x.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------
rng(trial);
n = length(x);
noise = (randn(n,1)+1j*randn(n,1))/sqrt(2);
sigma = (norm(x(:),'fro')/norm(noise(:),'fro'))*10^(-SNR/20);
sigma_2 = sigma^2;
z = sigma*noise;
y = x + z;
end