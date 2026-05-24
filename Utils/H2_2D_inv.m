function [X] = H2_2D_inv(Y)
% -------------------------------------------------------------------------
% This function is used for the pseudo-inverse 2D Hankel lifting operation.
% version 1.0 - 05/23/2026 (In this version, X must be a square matrix.)
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------------------------
N = sqrt(size(Y,1))*2-1; % N must be an odd number
N2 = (N+1)/2;
[W,~] = construct_W(N);

slice = zeros(N2,N2,2*N2-1);
indexs = zeros(1,N2*N2);
for i = 1:N2
    for j = 1:N2
        temp = i+j-1;
        slice(:,:,temp) = slice(:,:,temp) + Y((i-1)*N2+1:i*N2,(j-1)*N2+1:j*N2);
        indexs((i-1)*N2+j) = i+j-1;
    end
end

X = zeros(N,N);
for i = 1:2*N2-1
   slice(:,:,i) = slice(:,:,i) /sum(indexs==i);
   X(:,i) = H2_star(slice(:,:,i)).*W;
end

end