function [Y] = H2(x)
% -------------------------------------------------------
% This function is used for the Hankel lifting operation.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------
n = length(x);
n1 = (n + mod(n,2))/2;
Y = hankel(x(1:n1,1),x(n1:n,1));

end