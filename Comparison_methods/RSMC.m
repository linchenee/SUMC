function [x] = RSMC(y,Hx,n,omega,max_iter,rho,mu,tol,r,lambda,Idx)
% ----------------------------------------------------------    
% Reference-based structured matrix completion (RSMC) method
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ----------------------------------------------------------
[G,~,K] = svds(Hx,r);
prior = G*K';
max_mu = 1e7;
x = y;
L = zeros(size(Hx));
for iter = 1 : max_iter      
    % fprintf('%d\n',iter);
    [u,s,v] = svd(Hx-L/mu+lambda*prior/mu,'econ');
    N = u*max(s-1/mu,0)*v';

    x_last = x;
    x = H2_inv(N+L/mu,n,Idx);
    x(omega) = y(omega);
    Hx = H2(x);

    chg = norm(x - x_last,'fro')/norm(x_last,'fro');
    if chg < tol && iter >= 100 
       break;
    end 
    L = L + mu*(N-Hx);
    mu = min(rho*mu,max_mu);   
end

end