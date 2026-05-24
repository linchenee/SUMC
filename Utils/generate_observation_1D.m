function [xs, omega, x_star, f, amp] = generate_observation_1D(n, m, r, i)
% -----------------------------------------------------
% This function is used to generate the 1D observation.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -----------------------------------------------------
rng(round(i/2));
omega = randsample(n,m);
rng(i);
dynamic_range = 20;
c = exp(1i*2*pi*rand(r,1)).*10.^(rand(r,1)*dynamic_range/20);

xs = zeros(n,1);
freq_seed = randperm(1*n, r)/(1*n);
x_star = exp(2*pi*1i * (0:(n-1))' * freq_seed) * c;  
xs(omega) = x_star(omega);
[f, orders] = sort(angle(exp(1i*2*pi* freq_seed)));
f = f/(2*pi);
amp = c(orders).';

end