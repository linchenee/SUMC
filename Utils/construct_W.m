function [W,W_hat] = construct_W(n1)
% ------------------------------------------------------------------------------------------------
% This fucntion is used to construct the domain conversion matrices W and W_hat. 
% For details, see https://github.com/linchenee/ES-CPD/blob/main/functions/utilities/construct_W.m
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------------------------------------------------------------
if mod(n1,2) == 0
 W_hat = [1:n1/2,n1/2:-1:1]';
else
 W_hat = [1:(n1+1)/2,(n1-1)/2:-1:1]';
end
W = ones(n1,1)./W_hat;

end