function [x] = EMAC(y,n,omega,max_iter,rho,mu,tol,Idx)
% -----------------------------------------
% enhanced matrix completion (EMaC) method
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -----------------------------------------
max_mu = 1e7;
x = y;
Hx = H2(x);
A = randn(size(Hx));
for iter = 1 : max_iter
    [u,s,v] = svd(Hx-A/mu,'econ');
    N = u*max(s-1/mu,0)*v';

    x_last = x;
    x = H2_inv(N+A/mu,n,Idx);
    x(omega) = y(omega);
    Hx = H2(x);

    dY1 = N-Hx;     
    chg = norm(x - x_last,'fro')/norm(x_last,'fro');
    if chg < tol && iter >= 100
       break;
    end 
    A = A + mu*dY1;
    mu = min(rho*mu,max_mu);    
end

end