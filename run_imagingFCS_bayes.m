% This script loads the simulated fluorescence movie of a microdomain, runs 
% the automated blocking analysis, compute TACFs and covariance matrices, and
% runs the Bayesian analysis. This file was tested to run properly on MATLAB
% R2011a with Statistics Toolbox and Image Processing Toolbox installed.
% ----input parameter------------------------------------------------------
% o.model:  1D cell of model set to evaluate. See the function "FCS_fit_bayes_im"
%         for currently available models.
% o.sec_per_frame: sampling time of the movie
% o.corrlimit: maximal lag for TACFs
% o.um_per_px: pixel size in the image space, microns per pixel
% o.psf_sigma_um: standard deviation of the PSF in micron
% o.structure_f:  aspect ratio of the focal volume
% o.offset: camera offset (dark background)
% o.pbc: option for photobleaching correction (1--on, 0--off)
% o.plot:  option for plotting model fits (1--on, 0--off)
% o.binning: bin size for binning pixels in the movie (1 = no binning, 2 = 2X2 and so on) 
% o.fit:  option for nonlinear regression
%                       - 'GLS' (default): incorporate the noise correlations
%                       - 'w': weighted fitting (only noise level is used)
% o.regularization: Option for the shirnkage target of
% the regularized covariance matrix of the noise
%                       - '2way' (default): two-way target
%                       - 'none' : no regularization
%                       - 'D': Target D
%                       - 'B': Target B
% -----------------------------------------------------------------
% Copyright MIT 2013
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Oct 31, 2013
% -----------------------------------------------------------------
addpath('.\functions')
% load('.\microdomain_200k.mat');
%%
file_name  = 'C:\data\DiI longer acquisition time\simD2_seed1.tif' ;
im = readtiff(file_name) ;
%% inputs for the analysis
im = double(im) ; 
o.model = {'1comp2Ddiff_im' '2comp2Ddiff_im' '3comp2Ddiff_im'} ;  % fitting models 
o.sec_per_frame =  1e-3 ; 
% o.corrlimit =  5/o.sec_per_frame ; 
o.corrlimit =  Inf ; 
o.um_per_px = 0.4 ;
o.psf_sigma_um = 0.2 ;
o.offset = 0 ;
o.pbc = 0 ;
o.fit = 'GLS' ;
o.regularization = '2way' ;
% o.regularization = 'none' ;
o.binning = 1 ;
o.structure_f = 0 ; 
o.plot = 1 ;

% substract the camera offset
im = im - o.offset ;
% bin the pixels
if o.binning ==1
else
im = bin_image_3(im, o.binning) ;
end
imm = mean(im, 3);
%% photobleaching correction
if o.pbc ==1
    im = FCM_pb_correct(im,'pixelwise2') ;
else
end
%% run blocking analysis
[block_min_ind block_cur block_err t_block]= blocking_im(im, o.sec_per_frame, 1, o.corrlimit, 1) ;
    
%% use the maximal block-time for pixels failing the blocking test
block_min_ind_1 = block_min_ind ;
block_min_ind_1(block_min_ind == 0)=numel(t_block) ;

%% compute TACFs
j=1 ;
[Tcorr t cross_prod_binned]= compute_ACF_im(im,o.sec_per_frame, block_min_ind_1, o.corrlimit) ; % compute CCFs      
corrFCS_cell{j} = Tcorr;
cross_prod_cell{j} = cross_prod_binned ;        
%% plot computed TACFs
corrFCS = corrFCS_cell{1} ;
figure(9)    
for k = 1:size(corrFCS, 1)
plot(t, squeeze(corrFCS(k,:,:)),'linewidth',2)
hold all
end 
set(gca,'xscale','log');
xlabel('\tau (s)','FontSize',10)
ylabel('G (\tau)','FontSize',10)        
set(gca, 'xtick',10.^[-9:9])
axis tight
format_fig2(1)
%% perform Bayesian analysis
[fit_para noise_level]= imaging_FCS_bayes_analysis...
    (t, corrFCS_cell, cross_prod_cell, o) ;
%% generate maps
[mp D_map] = model_likelyhood_map(imm, fit_para, block_min_ind, o) ; 
