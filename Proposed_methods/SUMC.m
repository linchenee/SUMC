function [x] = SUMC(y,Hx,max_iter,n,omega,rho,mu_init,tol,r,w,ITERS_num,Idx)
% ------------------------------------------------    
% Subspace updated matrix completion (SUMC) method
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------------
max_mu = 1e7;
EYE = eye(size(Hx,1));
[n1,n2] = size(Hx);
x = y;

for ITERS = 1:ITERS_num
    [G, ~] = takagi(Hx,r);
    M = EYE - w*(G*G');
    L1 = zeros(n1,n2);
    L2 = zeros(n1,n2);
    L3 = zeros(n1,n2);
    x = y;
    D = Hx;
    B = M*D;
    temp1 = (conj(M)*M.'+EYE);
    temp2 = (M'*M+EYE);
    mu = mu_init;
    for iter = 1 : max_iter
        [u,s,v] = svd(B*conj(M)-L1/mu,'econ');
        A = u*max(s-1/mu,0)*v';
        B = ((A+L1/mu)*M.'+M*D-L2/mu)/temp1;
        D = temp2\(M'*(B+L2/mu)+Hx-L3/mu);

        x_last = x;
        x = H2_inv(D+L3/mu,n,Idx);
        x(omega) = y(omega);
        Hx = H2(x);
         
        chg = norm(x - x_last,'fro')/norm(x_last,'fro');
        if chg < tol && (iter >= 100 || ITERS>1)
           break;
        end 
    
        L1 = L1 + mu*(A-B*conj(M));
        L2 = L2 + mu*(B-M*D);    
        L3 = L3 + mu*(D-Hx);    
        mu = min(rho*mu,max_mu);    
    end
end

end