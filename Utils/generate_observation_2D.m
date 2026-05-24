function [xs,K,x_star] = generate_observation_2D(m,nd,r,c,delta_f,Ts,Doppl,Delay)
% -----------------------------------------------------
% This function is used to generate the 2D observation.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -----------------------------------------------------
N1 = nd(1);
N2 = nd(2);  
K = randsample(N1*N2,m);
X_star = zeros(N1,N2)+1j*zeros(N1,N2);
for n1 = 1:N1
    for n2 = 1:N2
        for l = 1:r
            X_star(n1,n2) = X_star(n1,n2)+c(l)*exp(-1i*2*pi*(n1-1)*delta_f*Delay(l)+1i*2*pi*(n2-1)*Ts*Doppl(l));
        end
    end
end
x_star = reshape(X_star,[N1*N2,1]);

xs = zeros(N1*N2,1);
xs(K) = x_star(K);
    
end

    
        
     