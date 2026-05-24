function [Hx] = init_2D(x,r)
% ----------------------------------------------------------------------------------------------
% This function is used to initialize the SUMC/SUMC-L method from the incomplete 2D observation.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ----------------------------------------------------------------------------------------------
[G, S] = takagi(H2_2D(x),r);
Hx = H2_2D(H2_2D_inv(G*S*G.'));

end