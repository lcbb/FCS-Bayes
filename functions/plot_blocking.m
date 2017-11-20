function plot_blocking(block_cur, block_err, t_block, block_min_ind)
% Plot blocking curves and its minimal block time

figure(84)   
for cur = 1:1         
%     for cur = 1:size(block_cur,1)         
        if block_min_ind(1,cur)>0
            errorbar(t_block, squeeze(block_cur(cur,:)),squeeze(block_err(cur,:)),...
                'linewidth',2, 'markersize',8);
            hold all
            
        else
            errorbar(t_block, squeeze(block_cur(cur,:)),squeeze(block_err(cur,:)),...
                'color', [0.7 0.7 0.7],'linewidth',2, 'markersize',8);    
            hold all
        end
    
    end
    axis tight
    xlim = get(gca, 'xlim') ; ylim = get(gca, 'ylim') ; 
    hold on
%     for cur = 1:size(block_cur,1)         
for cur = 1:1       
        if block_min_ind(cur)>0
            
            plot(t_block(block_min_ind(cur)), block_cur(cur,block_min_ind(cur)),...
                '.', 'markersize',40)
            hold all
%             area([0.9*t_block(block_min_ind_t(cur)) 1.1*t_block(block_min_ind_t(cur))], [ylim(2) ylim(2)]);    
%             hold all
        end
    
    end

    set(gca,'xscale','log');    
    set(gca, 'xtick',10.^[-10:10])  
    ylabel('\sigma','FontSize',10)
    xlabel('block time (s)','FontSize',10)
    format_fig2(5)
        hold off
%     print(gcf ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png'])
%     saveas(gcf,[fig_path, num2str(fig_num, '%03d'),'.fig'])
%      fig_num = fig_num +1 ;
end 
