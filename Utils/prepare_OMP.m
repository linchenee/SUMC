function [A,hatA] = prepare_OMP(omega,grid_num,n)
% ----------------------------------------------------   
% This function is used to prepare for the OMP method.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ----------------------------------------------------
Angle = (1:grid_num)/grid_num;
A = exp(2*pi*1i * (0:(n-1))' * Angle);
mask = zeros(n,1);
mask(omega) = 1;
mask2 = repmat(mask, [1,grid_num]);
hatA =  A.*mask2;

end