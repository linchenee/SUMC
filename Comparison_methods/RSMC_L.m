function [z] = RSMC_L(y,Hz,n,omega,max_iter,mu,tol,r,lambda,Idx)
% ------------------------------------------------    
% RSMC with low-rank factorization (RSMC-L) method
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------------
EYER = eye(r);

[G,S,K] = svds(Hz,r);
prior = G*K';
U = G*S;
V = K;
A = zeros(size(Hz));
z = y;
for iter = 1 : max_iter      
    z_last = z;
    z = H2_inv(U*V'-A,n,Idx);
    z(omega) = y(omega);
    Hz = H2(z);

    temp = (mu*(Hz+A)+lambda*prior);
    U = (temp*V)/(EYER+mu*(V'*V));
    V = (temp'*U)/(EYER+mu*(U'*U));
   
    chg = norm(z - z_last,'fro')/norm(z_last,'fro');
    if chg < tol && iter >= 100
        break;
    end 
    A = Hz-U*V' + A;
end
end