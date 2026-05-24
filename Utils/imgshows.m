function [] = imgshows(i,thr1,thr2,x,titles)
% -----------------------------------------
% The function is used to show results.
% version 1.0 - 05/23/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -----------------------------------------
figure(i);
set(gcf, 'Position', [100, 100, 900, 415]);
subfig = tight_subplot(2,4,[.145 .055],[.09 .04],[.05 .06]);      
for j = 1:8
    axes(subfig(j));
    imagesc(mean(squeeze(x(:,:,:,j)),3)); colormap(gray);%colormap(flipud(gray));
    set(gca,'YDir','normal','fontsize',9.5,'FontName','Times new roman');
    clim([thr1,thr2]);
    ylabel('Sample ratio ({\itp})','fontsize',9.5,'FontName','Times new roman');
    yticks([1, 6, 11, 16]);
    yticklabels({'0.1','0.3','0.5','0.7'});
    xlabel('Rank ({\itr})','fontsize',9.5,'FontName','Times new roman');
    xticks([1 5 10 15]);
    xticklabels({'2','10','20','30'});
    title(titles{1,j},'fontsize',10.5,'FontName','Times new roman');
end
colorbar('Position', [0.956 0.089 0.015 0.8736],...
    'Ticks',[0 0.5 1], 'FontSize',10, 'FontName','Times New Roman');
end