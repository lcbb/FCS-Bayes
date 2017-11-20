% This function runs the Bayesian analysis on FCS data with TACF curves and 
% their underlying independent intensity products.
%  - estimate noise and noise correlations from intensity products 
%  - calculate model probabilites and parameter estimates for each TACF curve 
%
% ----input------------------------------------------------------
% t:  column vector of the lag time
% corrFCS_cell:  1D cell of 3D arrays of TACF curves (pixel, pixel, lag time)
% cross_prod:  2D cell array with each cell recording the matrix of block-averaged intensity products for one pixel
% o.model:  1D cell of the model set to evaluate. See the function "FCS_fit_bayes_im"
%         for currently available models.
% o.sec_per_frame: sampling time of the movie
% o.corrlimit: maximal lag for TACFs
% o.um_per_px: pixel size in the image space, microns per pixel
% o.psf_sigma_um: standard deviation of the PSF in micron
% o.structure_f:  aspect ratio of the focal volume. 
% o.plot:  Option for plotting model fits (1--on, 0--off)
% o.binning: bin size for binning pixels in the movie (1 = no binning, 2 = 2X2 and so on) 
% o.fit:  Option for nonlinear regression
%                       - 'GLS' (default): incorporate the noise correlations
%                       - 'w': weighted fitting (only noise level is used)
% o.regularization: Option for the shirnkage target of
% the regularized covariance matrix of the noise
%                       - '2way' (default): two-way target
%                       - 'none' : no regularization
%                       - 'D': Target D
%                       - 'B': Target B
% ----Output-------------------------------------------------------
% fit_para: 2D, m-by-6 cell array of fitted parameters of m model fits
%           - column 1: parameter estimates of model fits (see model
%           functions for details)
%           - column 2: fitting uncertainties of parameters
%           - column 3: log model probabilities
%           - column 4: residuals of fits
%           - column 5: reduced chi-square values
%           - column 6: corresponding p-values of the chi-square value
% noise_level: 2D array of the noise level map
%-----------------------------------------------------------------
% Ref: 1. He, J. et al., Anal Chem 2012 
%      2. Guo, S. et al., Anal Chem 2012
%      3. Raftery, A. E., Sociological Methodology 1995
%      4. Sivia, D. S. et al., Data Analysis: A Bayesian Tutorial 2006
%-----------------------------------------------------------------
% Copyright MIT 2013
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Oct 31, 2012
% -----------------------------------------------------------------
%
function  [fit_para noise_level]= imaging_FCS_bayes_analysis(t, corrFCS_cell, cross_prod, o, varargin)
%% Check input options
% default options
model = o.model ;
ks = o.structure_f ;
plot_opt = o.plot ; 
fit_opt = o.fit ;
regularize_opt = o.regularization ;
n_model = size(model,2) ;
if nargin < 5, fix_opt = cell(n_model,1); fix_val = cell(n_model,1) ; end

curve_arr = corrFCS_cell{1} ;
prod_arr = cross_prod{1};
noise_level = zeros(size(curve_arr(:,:,1))) ;

clear cross_prod
n_cur = size(curve_arr,2) ;

fit_para = cell(n_model,6 , n_cur, size(curve_arr,1)) ;

h = waitbar(0,'Running Bayesian analysis...'); 
    for j=1:size(curve_arr,1)
%         close all
        corr = squeeze(curve_arr(j,:,:)) ;  
%% Plot input curves
        if plot_opt == 1
            figure(56)
            plot(t, corr)
            set(gca,'xscale','log');
            xlabel('\tau (s)','FontSize',10)
            ylabel('G (\tau)','FontSize',10)
            set(gca, 'xtick',10.^[-7:3])
            axis tight
            format_fig2(1)
        end

        fit_limit = floor(1*numel(t));
        corrFCS_1d = corr ;


    %% Noise level estimation
    cov_mat = zeros(numel(t), numel(t), n_cur) ; 
    ste_it  = zeros(size(curve_arr)) ;


    %% Estimate the noise covariance
    clear cross_prod_bin
    if strcmp(fit_opt, 'GLS')
        if strcmp(regularize_opt, 'none')
            for k = 1:n_cur
                cross_prod = prod_arr{j,k} ;          
                cov_mat(:,:,k) = cov(cross_prod)/size(cross_prod,1)  ;
                ste_it(j,k,:) = sqrt(diag(cov_mat(:,:,k))) ;
            end                 
%                 cov_mat = cov_it ;  
        elseif strcmp(regularize_opt, '2way')
%             cov_mat = zeros(size(cov_it)) ;            
            for k = 1:n_cur
                cross_prod = prod_arr{j,k} ;  
                [cov_mat(:,:,k), lamcor, lamvar] = covshrinkKPM(cross_prod, 1);
                cov_mat(:,:,k) = cov_mat(:,:,k)/size(cross_prod,1) ;
                ste_it(j,k,:) = sqrt(diag(cov_mat(:,:,k))) ;
            end
        else 
%             cov_mat = zeros(size(cov_it)) ;            

            for k = 1:n_cur
            cross_prod = prod_arr{j,k} ; 
            corr_dev =bsxfun(@minus,cross_prod, mean(cross_prod,1)) ;
            cov_mat(:,:,k) = corr_repair3(corr_dev, regularize_opt)/size(cross_prod,1);   
                ste_it(j,k,:) = sqrt(diag(cov_mat(:,:,k))) ;
            end            
        end
        
        if plot_opt == 1
            std_mat = sqrt(diag(cov_mat(:,:,1))) ;
            cov_mat_n = cov_mat(:,:,1)./(std_mat* std_mat') ;
            figure(10+j)
            imshow(cov_mat(:,:,1), 'InitialMagnification', 'fit')
            axis on
            colormap(jet)
            hc = colorbar ;
            set(gca, 'fontsize',15)
            caxis auto

            figure(20+j)
            imshow(cov_mat_n, 'InitialMagnification', 'fit')
            axis on
            colormap(jet)
            hc = colorbar ;
            set(gca, 'fontsize',15)
            caxis([-1 1])
          end

    sigma_i = permute(ste_it(j,:,1:fit_limit), [3 2 1]) ;

    end

    %% relative average noise level
    noise2sig = sigma_i(1:16,:)'./abs(corrFCS_1d(:,1:16)) ;
    noise_level(j,:) = mean(noise2sig,2) ;

    %% Fit individual curves and compute model probabilities
    corrFCS_good = corrFCS_1d ;
        for k=1:size(corrFCS_good,1)
        corrFCS_s = corrFCS_good(k,:)';
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
                        [a_1com std_beta logML resid_1com] ...
                        = FCS_fit_bayes_im(t(1:fit_limit),corrFCS_s(1:fit_limit),ks,'GLS',model{m},cov_mat(:,:,k), o.um_per_px*o.binning,o.psf_sigma_um(1),fix_opt{m} ,fix_val{m});

                        [chisq_1com pval_1com] = chisq_test(resid_1com, sigma_i(:,k), numel(a_1com)) ;
                        fit_para(m,:,k,j) = {a_1com std_beta logML resid_1com chisq_1com pval_1com} ;
                     end

                otherwise
                error('please specify the fitting option')
            end
            waitbar(j/size(curve_arr,1),h)
        %% Plot model fits and residuals
            if plot_opt == 1
                figure(100*j+k)
                h1 = subplot(2,1,1);    
                errorbar(t,corrFCS_s,sigma_i(:,k),...
                     'bo','markersize',2,'MarkerFaceColor','b','linewidth',1)
                t_extra = t;
                hold all

                for m = 1:size(model,2)
                    FCS_model_plot(fit_para{m,1,k,j},t_extra,ks,model{m},o.um_per_px*o.binning,o.psf_sigma_um(1))
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
                    plot(t ,fit_para{m,4,k,j}./sigma_i(:,k), 'linewidth',2)
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
    close(h)
end