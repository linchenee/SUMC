function [RMSE_ang,RMSE_amp,NRE] = evaluation(x0,x_star,freq,amp,n,r)
% --------------------------------------------------
% This function is used to evaluate the performance.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% --------------------------------------------------
temp = H2(x0);
f0 = rootmusic(temp*temp',r,'corr');
f0 = f0.'/(2*pi);
RMSE_ang = norm(sort(freq/pi*180)-sort(f0/pi*180),'fro')/sqrt(r);

amp0 = (pinv(exp(1i*2*pi*(0:(n-1))'*f0))*x0).';
RMSE_amp = norm(sort(amp0)-sort(amp),'fro')/sqrt(r);

NRE = norm(x0-x_star,'fro')/norm(x_star,'fro');

end