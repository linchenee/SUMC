function [Y] = H2_2D(X)
% ----------------------------------------------------------------------
% This function is used for the 2D Hankel lifting operation.
% version 1.0 - 05/23/2026 (In this version, X must be a square matrix.)
% Written by Lin Chen (lchen53@stevens.edu)
% ----------------------------------------------------------------------
N = size(X,1); % N must be a odd number
N2 = (N+1)/2;
for n = 1:N
   slice(:,:,n) = H2(X(:,n));
end
for i = 1:N2
    for j = 1:N2
        index = i+j-1;
        Y((i-1)*N2+1:i*N2,(j-1)*N2+1:j*N2) = slice(:,:,index);
    end
end

end