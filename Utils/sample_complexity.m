function [sam_EMaC,sam_RSMC,sam_SUMC,sam_log3,sam_log] = sample_complexity(x_star,n1,G,Prior,r,r_cor,w1,w2,Ak)   
% -----------------------------------------------------------------------------------
% This function is used to calculate the sample complexities of EMaC, RSMC, and SUMC.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -----------------------------------------------------------------------------------

n2 = n1; 
n = n1+n2-1;

[U,~] = takagi(H2(x_star),r);
mu1 = zeros(1,n1);
for i = 1:n1
    mu1(1,i) = norm(U(i,:),'fro')^2*n1/r;
end
mu1 = max(mu1); 

sam_EMaC = r*log(n)^4;
sam_log3 = r*log(n)^3;
sam_log = r*log(n);

[U0,~,V0] = svd(H2(x_star),'econ');
U0 = U0(:,1:r);
V0 = V0(:,1:r);

theta = w1*Prior;
F0_temp = U0*V0' - theta;

Left = (eye(n1)-U0*U0');
Right = (eye(n2)-V0*V0');
F0 = F0_temp - Left*F0_temp*Right;
PTc = Left*theta*Right;
F0_Ainf = 0;
F0_A2 = 0;
for k = 1:n
    temp1 = abs(sum(sum(F0.*Ak(:,:,k))));
    temp2 = norm(Ak(:,:,k));
    F0_Ainf = max(F0_Ainf,temp1*temp2);
    F0_A2 = F0_A2 + (temp1*temp2)^2;
end
F0_A2 = sqrt(F0_A2);
delta = 4*(F0_Ainf+F0_A2)/(1-2*norm(PTc));
sam_RSMC = max(delta^2,1) *r*log(n)^3*max(log(7*n*norm(F0,'fro')),1);

r_err = r - r_cor;
[~,~,R0] = svd(U'*G,'econ');
G_err = G*R0(:,r_cor+1:r);  
[Q1,R1] = qr(U*U'*G_err,0);
[Q2,R2] = qr((eye(n1) - U*U')*G_err,0);
[Q3,R3] = qr(eye(2*r_err) - w2*[R1;R2]*[R1;R2]');
for i = 1:2*r_err
 if R3(i,i) < 0
   R3(i,:) = -R3(i,:);
   Q3(i,:) = -Q3(i,:);
 end
end

mu2 = zeros(1,n1);
for i = 1:n1
 mu2(1,i) = norm(Q2(i,:),'fro')^2*n1/r;
end
mu2 = max(mu2);

P1 = R3(1:r_err,1:r_err);
P2 = R3(1:r_err,r_err+1:2*r_err);
delta2 = max(mu2/mu1,1)*(norm(P1)^4 + 2*norm(P1)^2*norm(P2)^2);
alpha1 = (7*n)^2/(exp(1)^3)*(r_cor*(1-w2)^4 + r_err*delta2);
alpha2 = (1-w2)^4 + delta2;
if alpha1 <= exp(1)
    f_delta2 = max(sqrt(alpha2)*log(n)^2,1);
else
    f_delta2 = log(n)^2*log(min(alpha1,n))*max(alpha2,1);
end
sam_SUMC = r*f_delta2*log(n);

end