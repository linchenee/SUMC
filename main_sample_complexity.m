%% Run this code to obtain the results in Fig. 2.
clc
clear
close all
addpath(genpath('.\Comparison_methods\'));
addpath(genpath('.\Proposed_methods\'));
addpath(genpath('.\Utils\'));

n1 = 50; % size of the Hankel matrix  
n2 = n1; % size of the Hankel matrix  
n = n1+n2-1; % length of the spectrally sparse signal 

r = 4; % rank 
r_cor = 2; % correct directions in the subspace knowledge
r_err = r - r_cor; % erroneous directions in the subspace knowledge
RMSEs_range = logspace(-1.25,-4.5,14); % tested RMSE
tol2 = 1e-10; % tolerance
w1 = 1; % weight hyperparameter for RSMC
w2 = 1; % weight hyperparameter for SUMC
Monte = 100; % Monte Carlo trials

Ak = zeros(n1,n2,n);
for k = 1:n
    ek = zeros(n,1);
    ek(k) = 1;
    temp = H2(ek);
    Ak(:,:,k) = temp/sqrt(sum(temp(:)));
end

sam_EMaC = zeros(length(RMSEs_range),Monte);
sam_RSMC = sam_EMaC;
sam_SUMC = sam_EMaC;
sam_log3 = sam_EMaC;
sam_log = sam_EMaC;
for i = 1:length(RMSEs_range)
    fprintf('%d\n',i);
    RMSE_i = RMSEs_range(i);
    for j = 1:Monte
        [~,~,x_star,freq,amp] = generate_observation_1D(n,1,r,j);
        mat_a = perms(freq.').';
        mat_amp = perms(amp.').';

        ids = randsample(r,r);
        f0_cor = freq(ids(1:r_cor));

        std_left = 1e-10; 
        std_right = 10; 
        noise = randn(1,r_err);
        for k = 1:2000
            std = (std_left + std_right)/2;
            f0_err = angle(exp(1i*2*pi*(freq(ids(r_cor+1:r))+std*noise)))/(2*pi);
            mat_b = repmat([f0_cor,f0_err].', 1, size(mat_a, 2));
            [RMSE1,order] = min(sqrt(sum(abs(mat_a/pi*180 - mat_b/pi*180).^2, 1)/r));
            f0 = mat_b(:,order).';
            f0_GT = mat_a(:,order).';
            amp_GT = mat_amp(:,order).';

            if abs(RMSE1 - RMSE_i) < tol2
                break;
            elseif RMSE1 < RMSE_i
                std_left = std;
            else 
                std_right = std;
            end
            if k == 2000
                fprintf('Divergent at i=%d, k=%d\n',i,j);
            end
        end
        
        Hx_est = H2(exp(2*pi*1i * (0:(n-1))' * f0) * amp_GT.'); 
        [U0, ~, V0] = svds(Hx_est,r);
        [G, ~] = takagi(Hx_est,r);
        Prior = U0*V0';
        [sam_EMaC(i,j),sam_RSMC(i,j),sam_SUMC(i,j),sam_log3(i,j),sam_log(i,j)] = sample_complexity(x_star,n1,G,Prior,r,r_cor,w1,w2,Ak);
    end
end

xlab = 1:length(RMSEs_range);
figure(1); % Fig. 2(a) in the paper
plot(xlab, mean(sam_EMaC,2),'g','LineWidth',3,'Marker','x','MarkerSize',9); hold on;
plot(xlab, mean(sam_RSMC,2),'b','LineWidth',3,'Marker','o','MarkerSize',9); hold on;
plot(xlab, mean(sam_SUMC,2),'r','LineWidth',3,'Marker','square','MarkerSize',9); hold on;
plot(xlab, mean(sam_log3,2),'k--','LineWidth',3,'MarkerSize',9); hold on;
plot(xlab, mean(sam_log,2),'k:','LineWidth',3,'MarkerSize',9); hold on;
xlabel('RMSE ( \circ )','fontsize',19,'FontName','Times new roman');
ylabel('Value','fontsize',19,'FontName','Times new roman');
xlim([1 14]);
ylim([0 2850]);
yticks([0 1000 2000]);
set(gca, 'FontName', 'Times new roman', 'FontSize', 19, 'TickLength',[0.02 0.025],...
    'XTick',[2 6 10 14],'XTickLabel',{'10^{-1.5}','10^{-2.5}','10^{-3.5}','10^{-4.5}'});
h = legend('$r \log^4(n)$ (EMaC)',...
    '$\delta_0 \cdot r \log^3(n)$ (RSMC)',...
    '$f(\delta_2) \cdot r \log(n)$ (SUMC)',...
    '$r \log^3(n)$','$r \log(n)$');
set(h,'Interpreter','latex','fontsize',16,'FontName','Times new roman',...
    'Position',[0.426398 0.672 0.459375 0.316286]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

n1 = 50; % size of the Hankel matrix  
n2 = n1; % size of the Hankel matrix  
n = n1+n2-1; % length of the spectrally sparse signal 

r = 8; % rank 
rs_cor = 0:7; % r_cor vector
rmse_ref = 10^(-3); % tested RMSE
tol = 1e-10; % tolerance
w1 = 1; % weight hyperparameter for RSMC
w2 = 1; % weight hyperparameter for SUMC
Monte = 100; % Monte Carlo trials

Ak = zeros(n1,n2,n);
for k = 1:n
    ek = zeros(n,1);
    ek(k) = 1;
    temp = H2(ek);
    Ak(:,:,k) = temp/sqrt(sum(temp(:)));
end

sam_EMaC = zeros(length(rs_cor),Monte);
sam_RSMC = sam_EMaC;
sam_SUMC = sam_EMaC;
sam_log3 = sam_EMaC;
sam_log = sam_EMaC;
for i = 1:length(rs_cor)
    fprintf('%d\n',i);
    r_cor = rs_cor(i);
    r_err = r - r_cor;
    for j = 1:Monte
        [~,~,x_star,freq,amp] = generate_observation_1D(n,1,r,j);
        mat_a = perms(freq.').';
        mat_amp = perms(amp.').';

        ids = randsample(r,r);
        f0_cor = freq(ids(1:r_cor));

        std_left = 1e-10; 
        std_right = 10; 
        noise = randn(1,r_err);
        for k = 1:2000
            std = (std_left + std_right)/2;
            f0_err = angle(exp(1i*2*pi*(freq(ids(r_cor+1:r))+std*noise)))/(2*pi);
            mat_b = repmat([f0_cor,f0_err].', 1, size(mat_a, 2));
            [rmse1,order] = min(sqrt(sum(abs(mat_a/pi*180 - mat_b/pi*180).^2, 1)/r));
            f0 = mat_b(:,order).';
            amp_GT = mat_amp(:,order).';

            if abs(rmse1 - rmse_ref) < tol
                break;
            elseif rmse1 < rmse_ref
                std_left = std;
            else 
                std_right = std;
            end
            if k == 2000
                fprintf('Divergent at i=%d, k=%d\n',i,j);
            end
        end
        
        Hx_est = H2(exp(2*pi*1i * (0:(n-1))' * f0) * amp_GT.'); 
        [U0, ~, V0] = svds(Hx_est,r);
        [G, ~] = takagi(Hx_est,r);
        Prior = U0*V0';
        [sam_EMaC(i,j),sam_RSMC(i,j),sam_SUMC(i,j),sam_log3(i,j),sam_log(i,j)] = sample_complexity(x_star,n1,G,Prior,r,r_cor,w1,w2,Ak);
    end
end

xlab = 1:length(rs_cor);
figure(2); % Fig. 2(b) in the paper
plot(xlab, mean(sam_EMaC,2),'g','LineWidth',3,'Marker','x','MarkerSize',9); hold on;
plot(xlab, mean(sam_RSMC,2),'b','LineWidth',3,'Marker','o','MarkerSize',9); hold on;
plot(xlab, mean(sam_SUMC,2),'r','LineWidth',3,'Marker','square','MarkerSize',9); hold on;
plot(xlab, mean(sam_log3,2),'k--','LineWidth',3,'MarkerSize',9); hold on;
plot(xlab, mean(sam_log,2),'k:','LineWidth',3,'MarkerSize',9); hold on;
xlabel('\itr_{\rm{cor}}','fontsize',19,'FontName','Times new roman');
ylabel('Value','fontsize',19,'FontName','Times new roman');
xlim([1 8]);
ylim([0 5650]);
set(gca,'FontName', 'Times new roman', 'FontSize', 19, 'TickLength',[0.02 0.025],...
    'XTick',[1 2 3 4 5 6 7 8],'XTickLabel',{'0','1','2','3','4','5','6','7'},...
    'YTick',[0 2000 4000]);
h = legend('$r \log^4(n)$ (EMaC)',...
    '$\delta_0 \cdot r \log^3(n)$ (RSMC)',...
    '$f(\delta_2) \cdot r \log(n)$ (SUMC)',...
    '$r \log^3(n)$','$r \log(n)$');
set(h,'Interpreter','latex','fontsize',16,'FontName','Times new roman',...
    'Position',[0.384 0.672 0.5 0.324]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

n1 = 50; % size of the Hankel matrix  
n2 = n1; % size of the Hankel matrix 
n = n1+n2-1; % length of the spectrally sparse signal

rs = 1:8; % rank vector
r_cor = 0; % correct directions in the subspace knowledge
rmse_ref = 1e-3; % tested RMSE
tol = 1e-10; % tolerance
w1 = 1; % weight hyperparameter for RSMC
w2 = w1; % weight hyperparameter for SUMC
Monte = 100; % Monte Carlo trials

Ak = zeros(n1,n2,n);
for k = 1:n
    ek = zeros(n,1);
    ek(k) = 1;
    temp = H2(ek);
    Ak(:,:,k) = temp/sqrt(sum(temp(:)));
end

sam_EMaC = zeros(length(rs),Monte);
sam_RSMC = sam_EMaC;
sam_SUMC = sam_EMaC;
sam_log3 = sam_EMaC;
sam_log = sam_EMaC;
for i = 1:length(rs)
    fprintf('%d\n',i);
    r = rs(i);
    r_err = r - r_cor;
    for j = 1:Monte
        [~,~,x_star,freq,amp] = generate_observation_1D(n,1,r,j);
        mat_a = perms(freq.').';
        mat_amp = perms(amp.').';

        ids = randsample(r,r);
        f0_cor = freq(ids(1:r_cor));

        std_left = 1e-10; 
        std_right = 10; 
        noise = randn(1,r_err);
        for k = 1:2000
            std = (std_left + std_right)/2;
            f0_err = angle(exp(1i*2*pi*(freq(ids(r_cor+1:r))+std*noise)))/(2*pi);
            mat_b = repmat([f0_cor,f0_err].', 1, size(mat_a, 2));
            [rmse1,order] = min(sqrt(sum(abs(mat_a/pi*180 - mat_b/pi*180).^2, 1)/r));
            f0 = mat_b(:,order).';
            amp_GT = mat_amp(:,order).';

            if abs(rmse1 - rmse_ref) < tol
                break;
            elseif rmse1 < rmse_ref
                std_left = std;
            else 
                std_right = std;
            end
            if k == 2000
                fprintf('Divergent at i=%d, k=%d\n',i,j);
            end
        end
        
        Hx_est = H2(exp(2*pi*1i * (0:(n-1))' * f0) * amp_GT.'); 
        [U0, ~, V0] = svds(Hx_est,r);
        [G, ~] = takagi(Hx_est,r);
        Prior = U0*V0';
        [sam_EMaC(i,j),sam_RSMC(i,j),sam_SUMC(i,j),sam_log3(i,j),sam_log(i,j)] = sample_complexity(x_star,n1,G,Prior,r,r_cor,w1,w2,Ak);
    end
end

xlab = 1:length(rs);
figure(3); % Fig. 2(c) in the paper
plot(xlab, mean(sam_EMaC,2),'g','LineWidth',3,'Marker','x','MarkerSize',9); hold on;
plot(xlab, mean(sam_RSMC,2),'b','LineWidth',3,'Marker','o','MarkerSize',9); hold on;
plot(xlab, mean(sam_SUMC,2),'r','LineWidth',3,'Marker','square','MarkerSize',9); hold on;
plot(xlab, mean(sam_log3,2),'k--','LineWidth',3,'MarkerSize',9); hold on;
plot(xlab, mean(sam_log,2),'k:','LineWidth',3,'MarkerSize',9); hold on;
xlabel('\itr','fontsize',19,'FontName','Times new roman');
ylabel('Value','fontsize',19,'FontName','Times new roman');
xlim([1 8]);
ylim([0 5500]);
set(gca,'FontName', 'Times new roman', 'FontSize', 19, 'TickLength',[0.02 0.025],...
    'XTick',[1 2 3 4 5 6 7 8],'XTickLabel',{'1','2','3','4','5','6','7','8'},...
    'YTick',[0 2000 4000]);
h = legend('$r \log^4(n)$ (EMaC)',...
    '$\delta_0 \cdot r \log^3(n)$ (RSMC)',...
    '$f(\delta_2) \cdot r \log(n)$ (SUMC)',...
    '$r \log^3(n)$','$r \log(n)$');
set(h,'Interpreter','latex','fontsize',16,'FontName','Times new roman',...
    'Position',[0.384 0.672 0.5 0.324]);