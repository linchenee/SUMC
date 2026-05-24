function [Hx] = init_1D(x,n,r,Idx)
% % -------------------------------------------------------------------------------------------------------    
% % This function is used to initialize the RSMC/SUMC/RSMC-L/SUMC-L method from the incomplete observation.
% % version 1.0 - 05/23/2026
% % Written by Lin Chen (lchen53@stevens.edu)
% % ------------------------------------------------------------------------------------------------------- 
[G, S] = takagi(H2(x),r);
Hx = H2(H2_inv(G*S*G.',n,Idx));

end