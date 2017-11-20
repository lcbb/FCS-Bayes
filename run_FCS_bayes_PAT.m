% This run file loads photon arrival time data, runs the automated blocking analysis, 
% compute TACFs and covariance matrices, and runs the Bayesian analysis. 
% This file was tested to run properly on MATLAB R2011a with Statistics 
% Toolbox and Image Processing Toolbox installed.
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Sep 13, 2012
% -----------------------------------------------------------------

addpath('.\functions')
model = {'1comp3Ddiff' 'trip1comp3Ddiff' 'trip2comp3Ddiff'};
dir_name ='D:\data\Thorsten_data\2011-08-01\PAT\Flu-LP\' ;
LP = {'1', '5', '10', '15', '20' , '30', '50'};
j = 6 ;
file_name = ['Flu-5nM-' LP{j} 'uW-'] ;

% options for "read_PAT" function. See the function for details
opt.bin = 2 ; % increase the binning size if you run into memory issues
opt.bits = 16 ;
opt.hdr = 1 ;
opt.n_trc = 2 ;

dt = 1/60e6* opt.bin ; % sampling time assuming the raw PAT has resolution of 1/60*10^6 sec

path_name = [dir_name file_name ,'1'] ;
[imA imB] = read_PAT(path_name, opt) ;
corrlimit = 0.1*size(imA,2);

ks = 6 ; % aspect ratio of the focal volume wz/wxy


% [block_cur block_err t_block] = blocking_analysis_ACF(imA, dt, 1, corrlimit) ; % blocking analysis for ACF
[block_cur block_err t_block] = blocking_analysis_CCF(imA,imB, dt, 1, corrlimit) ; % blocking analysis for CCF
block_min_ind = block_conver_cri41(block_cur, block_err) ; % find the minimal block time
plot_blocking(block_cur, block_err, t_block, block_min_ind) ; % plot the blocking curve and the minimal block time    
    
if block_min_ind == 0
    block_min_ind_s = numel(t_block) ;
else
    block_min_ind_s = block_min_ind ;
end

%     [Tcorr t cross_prod_binned]= compute_ACF2(imA, dt,  block_min_ind , corrlimit) ; % compute ACFs      
    [Tcorr t cross_prod_binned]= compute_CCF2(imA, imB, dt,  block_min_ind , corrlimit) ; % compute CCFs      
    corrFCS_cell{j} = Tcorr;
    cross_prod_cell{j} = cross_prod_binned ;     
    clear imA imB
    %% plot computed TACFs
    figure(9)
    plot(t, corrFCS_cell{j}, 'linewidth', 2)
    set(gca,'xscale','log');
    xlabel('\tau (s)','FontSize',10)
    ylabel('G (\tau)','FontSize',10)        
    set(gca, 'xtick',10.^[-9:9])
    axis tight
    format_fig2(1)

%% Compute model probabilities
[mp fit_para] = FCS_bayes_analysis_single_trace...
    (t, corrFCS_cell, cross_prod_cell, model, ks, 1, 'GLS', '2way') ;

[mp_m err_l err_h] = mp_stac(mp) ;  % get statistics of model probabilities

