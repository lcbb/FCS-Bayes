%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularized estimation of covariance matrix.
% A well-conditioned and more accurate estimator for covariance matrix when
% sample covariance is ill-conditioned (n~p)
% input: 2d array of observations with zero mean. Each row is one
% observation.
% output: 
% Ss - regularized estimated covariance matrix
% w1, w2 - coeffecients of the target matrix and the sample covariance matrix
% varargin - options for shrinkage targets
%           'I' - Indentity matrix defined in Ref 1.
%           'B' - Indentity matrix, target B defined in Ref 2.
%           'D' - Diagonal variance matrix, target D defined in Ref 2.
%
% Reference: 1. O. Ledoit and M. Wolf,J. Multivar. Anal., 2004 
%            2. J. Schäfer et. al., Statistical Applications in Genetics
%            and Molecular Biology, 2005
% 
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Aug.12.2011  
% -----------------------------------------------------------------

function [Ss w1 w2] = corr_repair3(corr_dev, varargin)
if isempty(varargin)
    option = 'D' ;
else 
    option = varargin{1} ;
end

p = size(corr_dev, 2) ; n = size(corr_dev, 1) ;

switch option
    case 'I'
        S = corr_dev'*corr_dev / n ;
        m =sum(diag(S))/p ;
        d2 = norm(S - m * eye(size(S)), 'fro')^2 /p ;

        % xtx = zeros(p,p) ;
        b2_bar = 0 ;
        for k = 1:n ;
        xtx = corr_dev(k,:)'*corr_dev(k,:) ;
        b2_bar = b2_bar + norm(xtx - S, 'fro')^2/p ;
        end

        b2_bar = b2_bar/n/n ; 
        b2 = min(b2_bar, d2) ;
        a2 = d2-b2 ;
        Ss = b2/d2 * m * eye(size(S)) + a2/d2 * S ;
        w1 = b2/d2 ; w2 = a2/d2 ;
        
    case 'B'
        w_bar = corr_dev'*corr_dev / n ;
        S = n/(n-1) * w_bar ;
        nu =sum(diag(S))/p ;
        var_s = 0 ;
        for k = 1:n ;
        wk = corr_dev(k,:)'*corr_dev(k,:) ;
        var_s = var_s + norm(wk - w_bar, 'fro')^2 ;
        end
        var_s = n/((n-1).^3)* var_s ;        
        lamda = var_s / norm(S - nu * eye(size(S)), 'fro')^2 ;
        w1 = lamda ; w2 = 1-w1 ;
        Ss = w1* nu * eye(size(S)) + w2*S ;
        
    case 'D'
        w_bar = corr_dev'*corr_dev / n ;
        S = n/(n-1) * w_bar ;
        var_s = 0 ;
        for k = 1:n ;
        wk = corr_dev(k,:)'*corr_dev(k,:) ;
        var_s = var_s + norm(wk - w_bar, 'fro')^2 - norm(diag(wk - w_bar),2)^2 ;
        end
        var_s = n/((n-1).^3)* var_s ;        
        lamda = var_s / (norm(S, 'fro')^2 - norm(diag(S),2)^2) ;
        w1 = lamda ; w2 = 1-w1 ;
        Ss = w1*diag(diag(S)) + w2*S ;
    otherwise
        error('Wrong option for regularization')
            
end

% [~, S_n] = cov2corr(S) ; [~, Ss_n] = cov2corr(Ss) ;

% figure(205)
% imshow(S_n, 'InitialMagnification', 'fit')
% axis on ; colormap(jet) ; h = colorbar ; set(gca, 'fontsize',15); caxis ([-1 1])
% 
% figure(206)
% imshow(Ss_n, 'InitialMagnification', 'fit')
% axis on ; colormap(jet) ; h = colorbar ; set(gca, 'fontsize',15); caxis ([-1 1])
% 
% figure(207)
% imshow(S, 'InitialMagnification', 'fit')
% axis on ; colormap(jet) ; h = colorbar ; set(gca, 'fontsize',15); caxis auto
% 
% figure(208)
% imshow(Ss, 'InitialMagnification', 'fit')
% axis on ; colormap(jet) ; h = colorbar ; set(gca, 'fontsize',15); caxis auto
end

