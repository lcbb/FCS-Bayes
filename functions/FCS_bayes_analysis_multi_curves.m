% This function runs the Bayesian analysis on FCS data with multiple TACF curves.
%  - estimate noise and noise correlations from multiple TACF curves 
%  - calculate model probabilites and parameter estimates for mean TACF curves 
% Note: All curves have to come from the same physical process.
%
% ----Required input------------------------------------------------------
% t:  column vector of the lag time
% corrFCS_1d:  2D array of TACF curves with individual curves in rows
% model:  1D cell of model set to evaluate. See the function "FCS_fit_bayes"
%         for currently available models.
% ----Optional input------------------------------------------------------
% ks:  aspect ratio of the focal volume. The default is 5.
% plot_opt:  Option for plotting model fits (1--on, 0--off)
% n_mean:  number of curves that each mean curve is calculated from. The 
%          default is calculating the mean curve from all input curves (generating only one mean curve)
% fit_opt:  Option for nonlinear regression
%                       - 'GLS' (default): incorporate the noise correlations
%                       - 'w': weighted fitting (only noise level is used)
% regularize_opt: Option for the shirnkage target of
% the regularized covariance matrix of the noise
%                       - '2way' (default): two-way target
%                       - 'none' : no regularization
%                       - 'D': Target D
%                       - 'B': Target B
% ----Output-------------------------------------------------------
% mp: 2D, k-by-m array of probabilities of m models for k mean curves
% fit_para: 2D, m-by-6 cell array of fitted parameters of m model fits
%           - column 1: parameter estimates of model fits (see model
%           functions for details)
%           - column 2: fitting uncertainties of parameters
%           - column 3: log model probabilities
%           - column 4: residuals of fits
%           - column 5: reduced chi-square values
%           - column 6: corresponding p-values of the chi-square value%
%-----------------------------------------------------------------
% Ref: 1. He, J. et al., Anal Chem 2012 
%      2. Guo, S. et al., Anal Chem 2012
%      3. Raftery, A. E., Sociological Methodology 1995
%      4. Sivia, D. S. et al., Data Analysis: A Bayesian Tutorial 2006
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 02, 2012
% -----------------------------------------------------------------
%
function [mp fit_para] = FCS_bayes_analysis_multi_curves(t, curve_cell, model, ks, plot_opt, n_mean, fit_opt, regularize_opt, fix_opt, fix_val)
%% Check input options
n_condi = size(curve_cell,2);
n_cur = size(curve_cell{1},1);
n_model = size(model,2) ;

% default options
if nargin < 10, fix_val = cell(n_model,1) ; end
if nargin < 9, fix_opt = cell(n_model,1) ; end
if nargin < 8, regularize_opt = '2way' ; end
if nargin < 7, fit_opt = 'GLS' ;  end
if nargin < 6, n_mean = n_cur ; end
if nargin < 5, plot_opt = 1 ;  end
if nargin < 4, ks = 5 ;  end

if mod(n_cur, n_mean) == 0  
else
   error('The number of mean curves must be a integer')
end

fit_para = cell(n_model,6 , n_cur/n_mean, n_condi);

for j=1:n_condi ;
    close all
    corrFCS_1d = curve_cell{j} ;
    %% Plot input curves
    if plot_opt == 1
        figure(56)
        plot(t, corrFCS_1d)
        set(gca,'xscale','log');
        xlabel('\tau (s)','FontSize',10)
        ylabel('G (\tau)','FontSize',10)
        set(gca, 'xtick',10.^[-7:3])
        axis tight
        format_fig2(1)
    end
    %% Noise level estimation

    corrFCS_std = repmat(std(corrFCS_1d),[n_cur 1]);
    if n_mean == 1
        corrFCS_mean = corrFCS_1d ;  
        corrFCS_ste = corrFCS_std ;    
    else    
        corrFCS_mean = binning_2d(corrFCS_1d',n_mean)'/n_mean ;
        corrFCS_ste = corrFCS_std/sqrt(n_mean);  
    end

    %% Estimate the noise covariance

    if strcmp(fit_opt, 'GLS')     
        if strcmp(regularize_opt, 'none')                                
               cov_mat = cov(corrFCS_1d)/n_mean;        
        elseif strcmp(regularize_opt, '2way')        
                cov_mat = covshrinkKPM(corrFCS_1d)/n_mean;              
        else                                  
                corr_dev =bsxfun(@minus,corrFCS_1d, mean(corrFCS_1d,1)) ;
                cov_mat = corr_repair3(corr_dev, regularize_opt)/n_mean;                  
        end
        corrFCS_ste = repmat(sqrt(diag(cov_mat))', [n_cur/n_mean 1]) ;    
        sigma_i = corrFCS_ste' ;
        cov_mat = repmat(cov_mat, [1 1 n_cur/n_mean]) ;  
    elseif strcmp(fit_opt, 'w')  
        sigma_i = corrFCS_ste' ;
    elseif strcmp(fit_opt, 'nw')  
        sigma_i = corrFCS_ste' ;
    else     
        error('Wrong fitting option. The fitting option must be "w" or "GLS"')
    end
    
%         [~, cov_mat_n] = cov2corr(cov_mat(:,:,1));
% 
% figure(12)
% imshow(cov_mat(:,:,1), 'InitialMagnification', 'fit')
%     axis on
%     colormap(jet)
%     h = colorbar ;
%     set(gca, 'fontsize',15)
%     caxis auto
% 
% figure(13)
% imshow(cov_mat_n, 'InitialMagnification', 'fit')
%     axis on
%     colormap(jet)
%     h = colorbar ;
%     set(gca, 'fontsize',15)
%     caxis([-1 1])
% tick = 10.^[-6:0] ;
% % tick = 10.^[-2:0] ;
% for k = 1:numel(tick)
% ind = find(t<=tick(k),1,'last') ;
% tick_itp(k) = ind+(tick(k)-t(ind))/(t(ind+1)-t(ind));
% end
% 
% tick_l = cellstr(num2str(log10(tick)','10^{%d}'));
% % tick_l = {'10_-6', '2', '4'} ;
% set(gca, 'xtick',tick_itp)
% set(gca, 'ytick',tick_itp)
% set(gca, 'XTickLabel',[])
% set(gca, 'yTickLabel',[])
% 
% xlTxt = text(tick_itp, 170*ones(size(tick_itp)), tick_l, ...   %# create text at same locations
%     'Interpreter','tex', ...                   %# specify tex interpreter
%     'VerticalAlignment','middle', ...             %# v-align to be underneath
%     'HorizontalAlignment','center',...
% 'fontsize',20);           %# h-aligh to be centered
% ylTxt = text(zeros(size(tick_itp))-1,tick_itp, tick_l, ...   %# create text at same locations
%     'Interpreter','tex', ...                   %# specify tex interpreter
%     'HorizontalAlignment','right',...
% 'fontsize',20);           %# h-aligh to be centered
% format_fig2(4)
% xlabel('\tau (s)','FontSize',30)
% ylabel('\tau (s)','FontSize',30)
    %% Fit mean curves and compute model probabilities
    for k=1:size(corrFCS_mean,1)
    corrFCS_s = corrFCS_mean(k,:)';
        switch fit_opt
            case 'w'
                for m = 1:size(model,2)
                    [a_1com std_beta logML resid_1com] ...
                    = FCS_fit_bayes(t,corrFCS_s,ks,'w', model{m},sigma_i(:,k),fix_opt{m}, fix_val{m});

                    [chisq_1com pval_1com] = chisq_test(resid_1com, sigma_i(:,k), numel(a_1com)) ;
                    fit_para(m,:,k,j) = {a_1com std_beta logML resid_1com chisq_1com pval_1com} ;
                end
                
            case 'nw'
                for m = 1:size(model,2)
                    [a_1com std_beta logML resid_1com] ...
                    = FCS_fit_bayes(t,corrFCS_s,ks,'nw', model{m},sigma_i(:,k),fix_opt{m}, fix_val{m});

                    [chisq_1com pval_1com] = chisq_test(resid_1com, sigma_i(:,k), numel(a_1com)) ;
                    fit_para(m,:,k,j) = {a_1com std_beta logML resid_1com chisq_1com pval_1com} ;
                end


           case 'GLS'
                 for m = 1:size(model,2)
                    [a_1com std_beta logML resid_1com ] ...
                    = FCS_fit_bayes(t,corrFCS_s,ks,'GLS',model{m},cov_mat(:,:,k),fix_opt{m}, fix_val{m});

                    [chisq_1com pval_1com] = chisq_test(resid_1com, sigma_i(:,k), numel(a_1com)) ;
                    fit_para(m,:,k,j) = {a_1com std_beta logML resid_1com chisq_1com pval_1com} ;
                 end

            otherwise
            error('Wrong fitting option. The fitting option must be "w" or "GLS"')
        end
    %% Plot model fits and residuals
        if plot_opt == 1
            figure(100+k)
            h1 = subplot(2,1,1);    
            errorbar(t,corrFCS_s,sigma_i(:,k),...
                 'bo','markersize',2,'MarkerFaceColor','b','linewidth',1)
            t_extra = t;
            hold all

            for m = 1:size(model,2)
                ph(m) = FCS_model_plot(fit_para{m,1,k,j},t_extra,ks,model{m}) ;
            end
            set(gca,'xscale','log');
            set(ph(size(model,2)), 'linestyle', '-.')
            xlabel('\tau (s)','FontSize',10)
            ylabel('G (\tau)','FontSize',10)    
            legend({'data' model{:}})
            set(gca, 'xtick',10.^[-7:3])
            axis tight
            format_fig2(7)
            hold off

            h2 = subplot(2,1,2);    
            hold all
            for m = 1:size(model,2)
                ph(m) = plot(t ,fit_para{m,4,k,j}./sigma_i(:,k), 'linewidth',2) ;
            end
            plot(t ,zeros(size(t)),'k:', 'linewidth',1)
            set(gca,'xscale','log');
            set(ph(size(model,2)), 'linestyle', '-.')
            xlabel('\tau (s)','FontSize',10)
            ylabel('resid/\sigma','FontSize',10)
            set(gca, 'xtick',10.^[-7:3])
            axis tight
            format_fig2(7)
            xlim = get(h1, 'xlim');
            set(gca, 'xlim', xlim)
            hold off       
        end
    end
end
%% Set NaN probablities to be zero  
for j = 1:n_condi          
    for k=1:n_cur/n_mean
             logML = cell2mat(fit_para(:,3,k,j));
             for g = 1:numel(logML)
                if isnan(logML(g))||~isreal(logML(g))
                    logML(g)= -Inf ;
                end
             end
            PM = exp(logML-max(logML)) ;
            PM = PM./sum(PM) ;        
            for m = 1:n_model
                mp{j}(k,m) = PM(m) ;        
            end
    end   
end
    %% Write model probabilities
%     dir_name = './' ;
%     file_name = 'MP' ;
%     fid=fopen([dir_name file_name '.txt'],'w');
%     model_name = char(model) ;
%     for m = 1:n_model
%         fprintf(fid,'%s\t', model_name(m,:));
%         for k = 1:n_cur/n_mean-1      
%             fprintf(fid,'%3.2f\t', mp(k,m));
%         end
%         fprintf(fid,'%3.2f\r\n', mp(end,m));     
%     end
%     fclose(fid);

end