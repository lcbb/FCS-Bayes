% This function fits FCS models to data and calculate log likelihoods of
% models
%
% ----Required input------------------------------------------------------
% time:  column vector of the lag time
% corr:  column vector of a single TACF curve 
% k: aspect ratio of the focal volume
% opt:  Option for nonlinear regression
%                       - 'GLS' : Generalized least squares
%                       - 'w': weighted least squares
%                       - 'nw': non-weighted least squares
% model:  Option for the model. Examples for model naming:
%         -'null': null model
%         -'1comp3Ddiff': 1-component 3D normal diffusion
%         -'1comp3Ddiff+flow': 1-component 3D normal diffusion plus a
%         uniform flow
%         -'2comp2D+3Ddiff': 2-component normal diffusion with one 2D
%         component and one 3D component
%         -trip1comp3Ddiff': 1-component 3D normal diffusion with triplet
%         blinking
%
%         Note: When adding new model to this function, modify the function
%         "FCS_model_plot" as well.
%
% noise:  column vector of the noise associated with the input TACF curve 
%
% ----Optional input(for fixing part of fitting parameters)-------------
% varargin{1}: binary row vector to specify the parameters to be fixed.
% varargin{2}: row vector to specify the values the parameters to be fixed at.
%              e.g. appending "[0 0 1 0 0], 5e-3" to the input for the model
%              '2comp3Ddiff' will fix the first diffusion time at 5 ms

% ----Output-------------------------------------------------------
% a:  vector of fitted parameters
% fit_para:  2D, m-by-6 cell array of fitted parameters of m model fits
% std_beta:  fitting uncertainties of parameters
% logML:  log model likelihood
% residual:  residuals of fits
% logCOVB:  contribution of parameter covariance to logML
% log_resid:  contribution of residuals to logML
% BIC:  BIC score of the model fit
%-----------------------------------------------------------------
% Ref: 1. He, J. et al., Anal Chem 2012 
%      2. Guo, S. et al., Anal Chem 2012
%      3. Raftery, A. E., Sociological Methodology 1995
%      4. Sivia, D. S. et al., Data Analysis: A Bayesian Tutorial 2006
% --------------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Jan 20, 2012
% -----------------------------------------------------------------
function [a std_beta logML residual logCOVB log_resid] = FCS_fit_bayes(time,corr,k,opt,model,noise, varargin)
%% Option for models
switch model
    case 'null'  % null model
    a0 = zeros(1,1);
    a0(1) = 0;
    f1 = @null_model ;
    
    case '1comp3Ddiff' % 1-component 3D diffusion
    a0 = zeros(1,3);
    a0(1) = (max(corr) - min(corr));
    a0(3) = max([min(corr) 0]);      
    [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
    a0(2) = time(t_ind) ;
%     a0(2) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
           
    f1 = @(a,t)diff1com3D(a,t,k) ;
    
    case '2comp3Ddiff'  % 2-component 3D diffusion
    a0 = zeros(1,5);        
    a0(1) = (max(corr) - min(corr))/2;
    a0(2) = a0(1);
    a0(5) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((2*a0(1)-a0(5))/2+ a0(5)))) ;
    a0(3) = time(t_ind) ;
    a0(4) = a0(3) ;
    f1 = @(a,t)diff2com3D(a,t,k) ;

    case '3comp3Ddiff'  % 3-component 3D diffusion
    a0 = zeros(1,7);
    a0(1) = (max(corr) - min(corr))/3; a0(2) = a0(1); a0(3) = a0(1);
    a0(4) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(5) = a0(4) ; a0(6) = a0(4) ;
    a0(7) = max([min(corr) 0]);
    f1 = @(a,t)diff3com3D(a,t,k) ;
         
    case '4comp3Ddiff'  % 4-component 3D diffusion
    a0 = zeros(1,9);
    a0(1) = (max(corr) - min(corr))/4; a0(2) = a0(1); a0(3) = a0(1); a0(4) = a0(1);
    a0(5) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(6) = a0(5) ; a0(7) = a0(5) ; a0(8) = a0(5) ;
    a0(9) = max([min(corr) 0]);
    f1 = @(a,t)diff4com3D(a,t,k) ;
    
    case 'flow' % 1-component flow
    a0 = zeros(1,3);
    a0(1) = (max(corr) - min(corr));
    a0(2) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(3) = max([min(corr) 0]);
    f1 = @flow1com ;
    
    case '1comp3Ddiff+flow' % 1-component 3D diffusion + flow
    a0 = zeros(1,3);
    a0(1) = (max(corr) - min(corr));
    a0(2) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(3) = a0(2)   ;
    a0(4) = max([min(corr) 0]);
    f1 = @(a,t)diff_flow1com3D(a,t,k) ;   
    
    case '2comp3Ddiff+flow' % 2-component 3D diffusion + flow
    a0 = zeros(1,5);        
    a0(1) = (max(corr) - min(corr))/2;  a0(2) = a0(1);
    a0(3) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(4) = a0(3); a0(5) = a0(3); 
    a0(6) = max([min(corr) 0]);
    f1 = @(a,t)diff_flow2com3D(a,t,k) ;

    case '1comp2Ddiff'
    a0 = zeros(1,3);
    a0(1) = (max(corr) - min(corr));
    a0(2) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(3) = max([min(corr) 0]);
    f1 = @diff1com2D ;
    
    case '2comp2Ddiff'
    a0 = zeros(1,5);        
    a0(1) = (max(corr) - min(corr))/2;  a0(2) = a0(1);
    a0(3) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(4) = a0(3);
    a0(5) = max([min(corr) 0]);
    f1 = @diff2com2D ;

    case '3comp2Ddiff'
    a0 = zeros(1,7);
    a0(1) = (max(corr) - min(corr))/3; a0(2) = a0(1); a0(3) = a0(1);
    a0(4) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(5) = a0(4) ; a0(6) = a0(4) ;
    a0(7) = max([min(corr) 0]);
    f1 = @diff3com2D ;
    
    case '1comp2Ddiff+flow' % 1-component 2D diffusion + flow
    a0 = zeros(1,3);
    a0(1) = (max(corr) - min(corr));
    a0(2) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(3) = a0(2)   ;
    a0(4) = max([min(corr) 0]);
    f1 = @diff_flow1com2D ;
    
    case '2comp2Ddiff+flow'
    a0 = zeros(1,5);        
    a0(1) = (max(corr) - min(corr))/2;  a0(2) = a0(1);
    a0(3) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(4) = a0(3); a0(5) = a0(3); 
    a0(6) = max([min(corr) 0]);
    f1 = @diff_flow2com2D ;
    
    case '2comp2D+3Ddiff' % 1-component 2D diffusion + 1-component 3D diffusion 
    a0 = zeros(1,5);        
    a0(1) = (max(corr) - min(corr))/2;  a0(2) = a0(1);
    a0(3) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(4) = a0(3);
    a0(5) = max([min(corr) 0]);
    f1 = @(a,t)diff2com2D_3D(a,t,k) ;
    
    case '3comp2_2D+1_3Ddiff'   % 2-component 2D diffusion + 1-component 3D diffusion
    a0 = zeros(1,7);
    a0(1) = (max(corr) - min(corr))/3; a0(2) = a0(1); a0(3) = a0(1);
    a0(4) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(5) = a0(4) ; a0(6) = a0(4) ;
    a0(7) = max([min(corr) 0]);
    f1 = @(a,t)diff3com_2_2D_1_3D(a,t,k) ;
    
    case '3comp1_2D+2_3Ddiff'   % 1-component 2D diffusion + 2-component 3D diffusion
    a0 = zeros(1,7);
    a0(1) = (max(corr) - min(corr))/3; a0(2) = a0(1); a0(3) = a0(1);
    a0(4) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(5) = a0(4) ; a0(6) = a0(4) ;
    a0(7) = max([min(corr) 0]);
    f1 = @(a,t)diff3com_1_2D_2_3D(a,t,k) ;
        
    case 'trip1comp3Ddiff' % 1-component 3D diffusion + triplet blinking
    a0 = zeros(1,5);
    
    % estimate the initial condition for high triplet fraction
%     n_bin = 20 ;
%     [nh,xout]= hist(corr,n_bin) ;
%     [nh_s ind] = sort(nh, 'descend') ;
%     a0(1) = xout(ind(2));
%     a0(3) = max([min(corr) 0]);
%     [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
%     a0(2) = time(t_ind) ;
%     a0(4) = 1-a0(1)/mean(corr(1:5)); 
    
       % estimate the initial condition for low triplet fraction
    a0(1) = (max(corr) - min(corr));
    a0(3) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
    a0(2) = time(t_ind) ;
    a0(4) = 0.15 ;
    a0(5) = 1e-6;
    f1 = @(a,t)trip_diff1com3D(a,t,k) ;
    
    case 'trip2comp3Ddiff'  % 2-component 3D diffusion + triplet blinking
    a0 = zeros(1,7);
    
    % estimate the initial condition for high triplet fraction
%     n_bin = 20 ;
%     [nh,xout]= hist(corr,n_bin) ;
%     [nh_s ind] = sort(nh, 'descend') ;
%     a0(1) = xout(ind(2))/2;
%     a0(2) = a0(1);
%     a0(5) = max([min(corr) 0]);
%     [dd t_ind] = min(abs(corr-((2*a0(1)-a0(3))/2+ a0(3)))) ;
%     a0(3) = time(t_ind) ;
%     a0(4) = a0(3) ;
%     a0(6) = 1-2*a0(1)/mean(corr(1:5)); 
    
    a0(1) = (max(corr) - min(corr))/2;
    a0(2) = a0(1);
    a0(5) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((2*a0(1)-a0(5))/2+ a0(5)))) ;
    a0(3) = time(t_ind) ;
    a0(4) = a0(3) ;
% a0(4) = 1e-3 ;
    a0(6) = 0.15 ;
    a0(7) = 1e-6;  
    f1 = @(a,t)trip_diff2com3D(a,t,k) ;
    
    case 'trip3comp3Ddiff'  % 3-component 3D diffusion + triplet blinking
    % estimate the initial condition for high triplet fraction
    a0 = zeros(1,9);
%     n_bin = 20 ;
%     [nh,xout]= hist(corr,n_bin) ;
%     [nh_s ind] = sort(nh, 'descend') ;
%     a0(1) = xout(ind(2))/3;
%     a0(2) = a0(1); a0(3) = a0(1);
%     a0(7) = max([min(corr) 0]);
%     [dd t_ind] = min(abs(corr-((3*a0(1)-a0(7))/2+ a0(7)))) ;
%     a0(4) = time(t_ind) ;
%     a0(5) = a0(4) ; a0(6) = a0(4) ;
%     a0(8) = 1-3*a0(1)/mean(corr(1:5));    

%     % estimate the initial condition for low triplet fraction        
    a0(1) = (max(corr) - min(corr))/3 ;
    a0(2) = a0(1); a0(3) = a0(1);
    a0(7) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((3*a0(1)-a0(7))/2+ a0(7)))) ;
    a0(4) = time(t_ind) ;
    a0(5) = a0(4) ; a0(6) = a0(4) ;
    a0(8) = 0.15;
    a0(9) = 1e-6;
    f1 = @(a,t)trip_diff3com3D(a,t,k) ;
    
    case 'trip4comp3Ddiff'  % 4-component 3D diffusion + triplet blinking
    a0 = zeros(1,11);
%     n_bin = 20 ;
%     [nh,xout]= hist(corr,n_bin) ;
%     [nh_s ind] = sort(nh, 'descend') ;
%     a0(1) = xout(ind(2))/4;
%     a0(2) = a0(1); a0(3) = a0(1); a0(4) = a0(1);
%     a0(9) = max([min(corr) 0]);
%     [dd t_ind] = min(abs(corr-((4*a0(1)-a0(9))/2+ a0(9)))) ;
%     a0(5) = time(t_ind) ;
%     a0(6) = a0(5) ; a0(7) = a0(5) ; a0(8) = a0(5) ;
%     a0(10) = 1-4*a0(1)/mean(corr(1:5));   

%     % estimate the initial condition for low triplet fraction        
    a0(1) = (max(corr) - min(corr))/4;
    a0(2) = a0(1); a0(3) = a0(1); a0(4) = a0(1);
    a0(9) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((4*a0(1)-a0(9))/2+ a0(9)))) ;
    a0(5) = time(t_ind) ;
    a0(6) = a0(5) ; a0(7) = a0(5) ; a0(8) = a0(5) ;
    a0(10) = 0.15;
    a0(11) = 1e-6;
    f1 = @(a,t)trip_diff4com3D(a,t,k) ; 
    
    case '2trip1comp3Ddiff' % 1-component 3D diffusion + triplet blinking
    a0 = zeros(1,7);
    
    % estimate the initial condition for high triplet fraction
%     n_bin = 20 ;
%     [nh,xout]= hist(corr,n_bin) ;
%     [nh_s ind] = sort(nh, 'descend') ;
%     a0(1) = xout(ind(2));
%     a0(3) = max([min(corr) 0]);
%     [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
%     a0(2) = time(t_ind) ;
%     a0(4) = 1-a0(1)/mean(corr(1:5)); 
    
       % estimate the initial condition for low triplet fraction
    a0(1) = (max(corr) - min(corr));
    a0(3) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
    a0(2) = time(t_ind) ;
    a0(4) = 0.15 ; a0(5)= a0(4);
    a0(6) = 1e-6;  a0(7)= a0(6);
    
    f1 = @(a,t)trip2_diff1com3D(a,t,k) ;
    
    case 'trip1comp3Ddiff_fk' % 1-component 3D diffusion + triplet blinking
    a0 = zeros(1,6);
    
    % estimate the initial condition for high triplet fraction
%     n_bin = 20 ;
%     [nh,xout]= hist(corr,n_bin) ;
%     [nh_s ind] = sort(nh, 'descend') ;
%     a0(1) = xout(ind(2));
%     a0(3) = max([min(corr) 0]);
%     [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
%     a0(2) = time(t_ind) ;
%     a0(4) = 5 ;  
%     a0(5) = 1-a0(1)/mean(corr(1:5));         
         % estimate the initial condition for low triplet fraction
    a0(1) = (max(corr) - min(corr));
    a0(3) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
    a0(2) = time(t_ind) ;
    a0(4) = 5 ;  
    a0(5) = 0.15 ;
    a0(6) = 1e-6;    
       
    f1 = @trip_diff1com3D_free_k ;
    
    case '1comp3Ddiff_fk'
    a0 = zeros(1,3);
    a0(1) = (max(corr) - min(corr));
    a0(2) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(3) = max([min(corr) 0]);
    a0(4) = 5 ;
    f1 = @diff1com3D_free_k ;
    
    case '2comp3Ddiff_fk'
    a0 = zeros(1,5);        
    a0(1) = (max(corr) - min(corr))/2;  a0(2) = a0(1);
    a0(3) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(4) = a0(3);
    a0(5) = max([min(corr) 0]);
    a0(6) = 2 ;
    f1 = @diff2com3D_free_k ;
    
    case '3comp3Ddiff_fk'
    a0 = zeros(1,7);
    a0(1) = (max(corr) - min(corr))/3; a0(2) = a0(1); a0(3) = a0(1);
    a0(4) = time(find(ismember(abs((max(corr) - min(corr))/2 - corr),min(min(abs((max(corr) - min(corr))/2 - corr)))),1,'first'));
    a0(5) = a0(4) ; a0(6) = a0(4) ;
    a0(7) = max([min(corr) 0]);
    a0(8) = 2 ;
    f1 = @diff3com3D_free_k ;
    
    case 'trip1comp3Drotdiff' % 1-component 3D diffusion + triplet blinking
    a0 = zeros(1,5);
    
    % estimate the initial condition for high triplet fraction
%     n_bin = 20 ;
%     [nh,xout]= hist(corr,n_bin) ;
%     [nh_s ind] = sort(nh, 'descend') ;
%     a0(1) = xout(ind(2));
%     a0(3) = max([min(corr) 0]);
%     [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
%     a0(2) = time(t_ind) ;
%     a0(4) = 1-a0(1)/mean(corr(1:5)); 
    
       % estimate the initial condition for low triplet fraction
    a0(1) = (max(corr) - min(corr));
    a0(3) = max([min(corr) 0]);
    [dd t_ind] = min(abs(corr-((a0(1)-a0(3))/2+ a0(3)))) ;
    a0(2) = time(t_ind) ;
    a0(4) = 0.15 ;
    a0(5) = 1e-6;
    a0(6) = 0.1 ;
    a0(7) = 10.^-4 ;
    f1 = @(a,t)trip_rot_diff1com3D(a,t,k) ;
    
    
end 

%% Fix certain fitting parameters
    if ~isempty(varargin)
       fix_opt = varargin{1} ;
       fix_opt = logical(fix_opt) ;
       fix_val  = varargin{2} ;
       a0(fix_opt) = [] ;
    f = @(a,t)localfit(@(b,t1)f1(b,t1),a,t,fix_opt, fix_val) ;
    else
        f = f1 ;
    end

    if ~strcmp(opt, 'GLS')    
        w = 1./ noise.^2;
    else    
        covrpr = noise ;
    end
          boxsize = 200; % flat prior box size
    if strcmp(opt, 'w')
        [a, residual,J,COVB,mse,sigma_sqi] = nlinfitw(time,corr,f,a0,w); %[bFitw,rw,Jw,COVBw,msew]       
        log_resid = - 0.5*( sum(residual.^2.*w));
        log_DL = log_resid - numel(time)/2*log(2*pi)- sum(log(sqrt(w))) ;
        
    elseif strcmp(opt, 'nw')
        [a, residual,J,COVB,mse] = nlinfit(time,corr,f,a0);
        w = ones(size(time))*(1/mean(1./w));       
        log_resid = - 0.5*( sum(residual.^2.*w));
        log_DL = log_resid - numel(time)/2*log(2*pi)- sum(log(sqrt(w))) ;
        
    elseif strcmp(opt, 'GLS')
        [a, residual,J,COVB] = nlinfit_GLS(time,corr,f,a0, covrpr);        
        log_resid = - 0.5*( residual' / covrpr * residual);
        log_DL = log_resid - numel(time)/2*log(2*pi)- 0.5*log(det(covrpr)) ;    % data marginal likelihood
    end
    
        std_beta = sqrt((diag(COVB)))'; % uncertainties of fitting parameters
        logML = 0.5*numel(a0)*log(2*pi)+ 0.5*log(det(COVB)) + log_resid;  % Laplace approximation
        logML = logML - log(prod(std_beta)*(2*boxsize)^numel(a));   % model probability
        logCOVB = 0.5*log(det(COVB));
        
   

%% Fill in fixed parameters
    if ~isempty(varargin)
        b = zeros(size(fix_opt)) ;
        b_std = zeros(size(fix_opt)) ;
        b(fix_opt) = fix_val ;
        b(~fix_opt) = a ;
        b_std(~fix_opt) = std_beta;
        a = b ;
        std_beta = b_std;        
    else
    end
   
end 

function y=localfit(fun, beta,x, fixed, bfixed)
   
   b(fixed) = bfixed;
   b(~fixed) = beta;
   y = fun(b,x);
end