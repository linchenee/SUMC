function [result] = SUMC_2D(Y,Hx,max_iter,rho,mu_init,tol,r,w,ITERS_num,eps)
% -------------------------------------------------------------------    
% Subspace updated matrix completion (SUMC) method in the 2D scenario
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------------------
Obs_Hankel = H2_2D(Y);
omega = find(Obs_Hankel);
max_mu = 1e7;
EYE = eye(size(Hx,1));
[n1,n2] = size(Hx);

for ITERS = 1:ITERS_num
    [G, ~] = takagi(Hx,r);
    M = EYE - w*(G*G');
    L1 = zeros(n1,n2);
    L2 = zeros(n1,n2);
    L3 = zeros(n1,n2);
    D = Hx;
    B = M*D;
    temp1 = (conj(M)*M.'+EYE);
    temp2 = (M'*M+EYE);
    mu = mu_init;
    for iter = 1 : max_iter
        % fprintf('%d\n',iter);
        [u,s,v] = svd(B*conj(M)-L1/mu,'econ');
        A = u*max(s-1/mu,0)*v';
        B = ((A+L1/mu)*M.'+M*D-L2/mu)/temp1;
        D = temp2\(M'*(B+L2/mu)+Hx-L3/mu);
       
        Hx_last = Hx;
        Hx = D+L3/mu;
        if norm(Hx(omega)-Obs_Hankel(omega)) >= eps
          Hx(omega) = Obs_Hankel(omega) - eps*(Obs_Hankel(omega)-Hx(omega))/norm(Obs_Hankel(omega)-Hx(omega));
        end
        
        Hx = H2_2D(H2_2D_inv(Hx));

        chg = norm(Hx - Hx_last,'fro')/norm(Hx_last,'fro');
        if chg < tol
           break;
        end 
    
        L1 = L1 + mu*(A-B*conj(M));
        L2 = L2 + mu*(B-M*D);    
        L3 = L3 + mu*(D-Hx);    
        mu = min(rho*mu,max_mu);    
    end
end

result = H2_2D_inv(Hx);

end