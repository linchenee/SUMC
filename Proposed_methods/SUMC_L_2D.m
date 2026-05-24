function [result] = SUMC_L_2D(Y,Hx,max_iter,rho,mu_init,tol,r,w,ITERS_num,eps)
% -------------------------------------------------------------------    
% SUMC with low-rank factorization (SUMC-L) method in the 2D scenario
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------------------
Obs_Hankel = H2_2D(Y);
omega = find(Obs_Hankel);
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
    mu = mu_init;
    for iter = 1 : max_iter
        Z = (tempA+2*mu*EYE)\(mu*(A-L1/mu+B-L2/mu));
        A = ((Hx-L3/mu)*conj(B)+(Z+L1/mu))/(B.'*conj(B)+EYER);
        B = ((Hx-L3/mu).'*conj(A)+(Z+L2/mu))/(A.'*conj(A)+EYER);
        ABT = A*B.';
        
        Hx_last = Hx;
        Hx = ABT+L3/mu;
        if norm(Hx(omega)-Obs_Hankel(omega)) >= eps
          Hx(omega) = Obs_Hankel(omega) - eps*(Obs_Hankel(omega)-Hx(omega))/norm(Obs_Hankel(omega)-Hx(omega));
        end

        Hx = H2_2D(H2_2D_inv(Hx));
        chg = norm(Hx - Hx_last,'fro')/norm(Hx_last,'fro');
        if chg < tol
            break;
        end 
        L1 = L1 + mu*(Z-A);
        L2 = L2 + mu*(Z-B);
        L3 = L3 + mu*(ABT-Hx);
        mu = min(rho*mu,max_mu);    
    end
end

result = H2_2D_inv(Hx);

end