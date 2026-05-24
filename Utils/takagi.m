function [U,S] = takagi(A,q)
% [U,S]=takagi(A,q) returns q largest singular values in decreasing order
% and the corresponding vectors of the Takagi decomposition of A; all
% singular values and singular vectors are computed using svd, and only the
% first q are retained
[U,S,~] = svd(A,'vector');
E2 = diag(U'*A*conj(U));
E2 = E2./S;
E = sqrt(E2);
U = U*diag(E);
if nargin == 1
    S = diag(S);
else
    U = U(:,1:q);
    S = diag(S(1:q));
end