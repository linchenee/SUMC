function [x] = SUMC_L(y,Hx,max_iter,n,omega,rho,mu_init,tol,r,w,ITERS_num,Idx)
% ------------------------------------------------    
% SUMC with low-rank factorization (SUMC-L) method
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------------
max_mu = 1e7;
EYE = eye(size(Hx,1));
EYER = eye(r);
for ITERS = 1:ITERS_num
    [G, S] = takagi(Hx,r);
    A = G*sqrt(S);
    B = A;
    L1 = zeros(size(A));
    L2 = zeros(size(A));
    L3 = zeros(size(Hx));
    M = EYE - w*(G*G');
    tempA = 2*(M'*M);
    x = y;
    mu = mu_init;
    for iter = 1 : max_iter
        Z = (tempA+2*mu*EYE)\(mu*(A-L1/mu+B-L2/mu));
        A = ((Hx-L3/mu)*conj(B)+(Z+L1/mu))/(B.'*conj(B)+EYER);
        B = ((Hx-L3/mu).'*conj(A)+(Z+L2/mu))/(A.'*conj(A)+EYER);
        ABT = A*B.';
        
        x_last = x;
        x = H2_inv(ABT+L3/mu,n,Idx);
        x(omega) = y(omega);
        Hx = H2(x);

        chg = norm(x - x_last,'fro')/norm(x_last,'fro');
        if chg < tol && (iter >= 100 || ITERS>1)
            break;
        end 
        L1 = L1 + mu*(Z-A);
        L2 = L2 + mu*(Z-B);
        L3 = L3 + mu*(ABT-Hx);
        mu = min(rho*mu,max_mu);    
    end
end

end