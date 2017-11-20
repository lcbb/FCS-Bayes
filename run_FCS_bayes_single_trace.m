% This script loads simulated FCS photon-count trace data sets of two-component
% diffusion with vaying D2 for each data set, runs the blocking analysis,
% computes TACFs and covariance matrices and performs the Bayesian analysis. 
% Each data set has 8 photon-count traces.
% This file was tested to run properly on MATLAB
% R2011a with Statistics Toolbox and Image Processing Toolbox installed.
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 23, 2012
% -----------------------------------------------------------------

addpath('.\functions')
model = {'1comp3Ddiff' '2comp3Ddiff' '3comp3Ddiff'};


D2 = 10.^linspace(0.8,2.8,15) ; % condtions for D2
D2 = D2(end:-1:1) ;

% pre-allocate variables
block_min_ind = cell(numel(D2), 1) ; block_min_ind_s = cell(numel(D2), 1) ;
block_cur= cell(numel(D2), 1) ; 
block_err = cell(numel(D2), 1) ; 
t_block = cell(numel(D2), 1) ;
corrFCS_cell = cell(numel(D2), 1) ;
cross_prod_cell = cell(numel(D2), 1) ;

ks = 2.0575 ; % simulation parameter for wz/wxz

% It might take hours to run all the 15 conditions. Change the FOR loop if
% you want to run less condtions.
for j = 1: numel(D2) 
    load(['.\examples\FCS_8int_vary_Dr', '_D2=', num2str(D2(j), '%.3g'), '.mat']) ; 
    corrlimit = 0.1*size(im,2) ;
    
   %% blocking analysis
    [block_cur{j} block_err{j} t_block{j}] = blocking_analysis_ACF(im, dt, 1, corrlimit) ; % blocking analysis for ACF
    block_min_ind{j} = block_conver_cri41(block_cur{j}, block_err{j}) ;
    plot_blocking(block_cur{j}, block_err{j}, t_block{j}, block_min_ind{j}) ;
    %% use the maximal block-time for intensity traces failing the blocking test
    block_min_ind_s{j} = block_min_ind{j} ;
    for k = 1:numel(block_min_ind_s{j})
       if block_min_ind_s{j}(k) == 0
          block_min_ind_s{j}(k) = numel(t_block{j}) ;
       else
       end
    end 
    %% compute TACFs
    [Tcorr t cross_prod_binned]= compute_ACF2(im, dt, block_min_ind_s{j}, corrlimit) ; % compute TACFs 
    corrFCS_cell{j} = Tcorr;
    cross_prod_cell{j} = cross_prod_binned ;     
    clear im
    %% plot computed TACFs
    figure(9)
    plot(t, corrFCS_cell{j}, 'linewidth', 2)
    set(gca,'xscale','log');
    xlabel('\tau (s)','FontSize',10)
    ylabel('G (\tau)','FontSize',10)        
    set(gca, 'xtick',10.^[-9:3])
    axis tight
    format_fig2(1)
end
%% Compute model probabilities
[mp fit_para] = FCS_bayes_analysis_single_trace...
    (t, corrFCS_cell, cross_prod_cell, model, ks, 0, 'GLS', '2way') ;

[mp_m err_l err_h] = mp_stac(mp) ;  % get statistics of model probabilities

%% plot model probabilities
Dr = D2./63.1;  % simulation parameter for D-ratio
x1 = Dr ;
l_style = {'rs-','bo-', 'k^-', 'g*-', 'cv-'} ;
l_color = {'r','b', 'k', 'g', 'c'};
n_model = numel(model) ;
figure(11)
for k =1:n_model
errorbar(x1, mp_m(k,:), err_l(k,:),err_h(k,:),...
    l_style{k},'MarkerFaceColor',l_color{k},'markersize', 10,'linewidth', 2)
hold on
end
set(gca,'xscale','log');
tick = 10.^[-7:3] ;
set(gca, 'xtick',tick)
axis tight
set(gca, 'ylim', [0 1]) ;
legend('N_D=1','N_D=2','N_D=3')
legend('location', 'east')
xlabel('D_2/D_1')
ylabel('P(M)','FontSize',15)
format_fig2(5)
hold off
