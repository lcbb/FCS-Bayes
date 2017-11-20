% This funtion performs generalized least square nonlinear regression
% input:
% x,y,modelFun,beta0 - same as the input for nlinfit
% covrpr - covariance matrix of errors
% output - refer to the output of nlinfit
% When covariance matrix is not semi-difinite, the matrix needs to be
% "repaired" for Cholesky decomopsition
% 
% References:
%      [1] Seber, G.A.F, and Wild, C.J. (1989) Nonlinear Regression, Wiley.
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Oct 28, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bFitw,resid,Jw,COVBw] = nlinfit_GLS(x,y,modelFun,beta0,covrpr)

options = statset('MaxIter' , 2000, 'DerivStep',10^-8);
L = chol(covrpr, 'lower') ;
z = L\y;
modelFun_k = @(b,x) L \ modelFun(b,x);
[bFitw,rw,Jw,COVBw] = nlinfit(x,z,modelFun_k, beta0, options);

resid = y - modelFun(bFitw,x);