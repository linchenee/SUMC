function [x] = H2_inv(Hx,n,Idx)
% ----------------------------------------------------------------------
% This function is used for the pseudo-inverse Hankel lifting operation.
% version 1.0 - 05/23/2026 (In this version, X must be a square matrix.)
% Written by Lin Chen (lchen53@stevens.edu)
% ----------------------------------------------------------------------
[W,~] = construct_W(n);
x = accumarray(Idx,Hx(:)).*W;

end