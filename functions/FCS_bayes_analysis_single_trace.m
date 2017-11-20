% This function runs the Bayesian analysis on FCS data with TACF curves and 
% their underlying independent intensity products.
%  - estimate noise and noise correlations from intensity products 
%  - calculate model probabilites and parameter estimates for each TACF curve 
%
% ----Required input------------------------------------------------------
% t:  column vector of the lag time
% corrFCS_1d:  2D array of TACF curves with individual curves in rows
% cross_prod:  independent intensity products given by the block-transformation
% model:  1D cell of model set to evaluate. See the function "FCS_fit_bayes"
%         for currently available models.
% ----Optional input------------------------------------------------------
% ks:  aspect ratio of the focal volume. The default is 5.
% plot_opt:  Option for plotting model fits (1--on, 0--off)
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
%           - column 6: corresponding p-values of the chi-square value
%-----------------------------------------------------------------
% Ref: 1. He, J. et al., Anal Chem 2012 
%      2. Guo, S. et al., Anal Chem 2012
%      3. Raftery, A. E., Sociological Methodology 1995
%      4. Sivia, D. S. et al., Data Analysis: A Bayesian Tutorial 2006
%-----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Sep 11, 2013
% -----------------------------------------------------------------
%

function [mp fit_para] = FCS_bayes_analysis_single_trace(t, curve_cell, cross_prod_cell, model, ks, plot_opt, fit_opt, regularize_opt)
%% Check input options
n_condi = numel(curve_cell);
n_cur = size(curve_cell{1},1);
n_model = size(model,2) ;

% default options
if nargin < 8, regularize_opt = '2way' ; end
if nargin < 7, fit_opt = 'GLS' ;  end
if nargin < 6, plot_opt = 1 ;  end
if nargin < 5, ks = 5 ;  end

fit_para = cell(n_model,6 , n_cur, n_condi) ;

for j=1:n_condi ;
%     close all
    corrFCS_1d = curve_cell{j} ;
    prod_arr = cross_prod_cell{j} ;
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
    ste_it = zeros(size(corrFCS_1d)) ;
    cov_it = zeros(size(corrFCS_1d,2), size(corrFCS_1d,2), size(corrFCS_1d,1)) ;
    for k = 1:n_cur
        cross_prod = prod_arr{k} ;
        cov_it(:,:,k) = cov(cross_prod)/size(cross_prod,1) ;
        ste_it(k,:) = sqrt(diag(cov_it(:,:,k))) ;
    end

    %% Estimate the noise covariance
    if strcmp(fit_opt, 'GLS')
        if strcmp(regularize_opt, 'none')
            cov_mat = cov_it ;  
        elseif strcmp(regularize_opt, '2way')
            cov_mat = zeros(size(cov_it)) ;               
            for k = 1:n_cur
                cross_prod = prod_arr{k} ;
                [cov_mat(:,:,k), lamcor, lamvar] = covshrinkKPM(cross_prod, 1);
                cov_mat(:,:,k) = cov_mat(:,:,k)/size(cross_prod,1) ;
                ste_it(k,:) = sqrt(diag(cov_mat(:,:,k))) ;
            end
        else 
            cov_mat = zeros(size(cov_it)) ;                    
            for k = 1:n_cur
                cross_prod = prod_arr{k} ;
                corr_dev =bsxfun(@minus,cross_prod, mean(cross_prod,1)) ;
                cov_mat(:,:,k) = corr_repair3(corr_dev, regularize_opt)/size(cross_prod,1);    
                ste_it(k,:) = sqrt(diag(cov_mat(:,:,k))) ;
            end     
        end
        
        if plot_opt == 1
        [~ , cov_mat_n] = cov2corr(cov_mat(:,:,1));
        figure(10+j)
        imshow(cov_mat(:,:,1), 'InitialMagnification', 'fit')
        axis on
        colormap(jet)
        h = colorbar ;
        set(gca, 'fontsize',15)
        caxis auto

        figure(20+j)
        imshow(cov_mat_n, 'InitialMagnification', 'fit')
        axis on
        colormap(jet)
        h = colorbar ;
        set(gca, 'fontsize',15)
        caxis([-1 1])
        end
    end     
    sigma_i = ste_it' ;
    %%
   
    %% Fit individual curves and compute model probabilities
    for k=1:size(corrFCS_1d,1)
    corrFCS_s = corrFCS_1d(k,:)';
         switch fit_opt
            case 'w'
                for m = 1:size(model,2)
                    [a_1com std_beta logML resid_1com] ...
                    = FCS_fit_bayes(t,corrFCS_s,ks,'w', model{m},sigma_i(:,k));

                    [chisq_1com pval_1com] = chisq_test(resid_1com, sigma_i(:,k), numel(a_1com)) ;
                    fit_para(m,:,k,j) = {a_1com std_beta logML resid_1com chisq_1com pval_1com} ;
                end

           case 'GLS'
                 for m = 1:size(model,2)
                    [a_1com std_beta logML resid_1com ] ...
                    = FCS_fit_bayes(t,corrFCS_s,ks,'GLS',model{m},cov_mat(:,:,k));

                    [chisq_1com pval_1com] = chisq_test(resid_1com, sigma_i(:,k), numel(a_1com)) ;
                    fit_para(m,:,k,j) = {a_1com std_beta logML resid_1com chisq_1com pval_1com} ;
                 end

            otherwise
            error('Wrong fitting option. The fitting option must be "w" or "GLS"')
        end
    %% Plot model fits and residuals
        if plot_opt == 1
            figure(100*j+k)
            h1 = subplot(2,1,1);    
            errorbar(t,corrFCS_s,sigma_i(:,k),...
                 'bo','markersize',2,'MarkerFaceColor','b','linewidth',1)
            t_extra = t;
            hold all

            for m = 1:size(model,2)
                FCS_model_plot(fit_para{m,1,k,j},t_extra,ks,model{m})
            end
            set(gca,'xscale','log');
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
                plot(t ,fit_para{m,4,k,j}./sigma_i(:,k), 'linewidth',1)
            end
            plot(t ,zeros(size(t)),'k:', 'linewidth',1)
            set(gca,'xscale','log');
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
    for k=1:n_cur
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