function [Y] = H2_star(X)
% -------------------------------------------------------------------------------
% This function is used to perform the adjoint of the 2D Hankel lifting operation.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------------------------------
[n1,n2,n3] = size(X);
Y = zeros(n1+n2-1,n3);
for k = 1:n3
 for i = 0:n1-1
  for j = 0:n2-1
    Y(i+j+1,k) = Y(i+j+1,k) + X(i+1,j+1,k); 
  end
 end
end
end

