% calculate reduced chi-square and p-value
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function [chisq_res_r p_value] = chisq_test(resid, sigma_i, n_para )

resid_norm = (resid'./sigma_i');
chisq_res = sum(resid_norm.^2,2);
chisq_res_r = chisq_res/(size(resid,1)-n_para) ;  % degree of freedom = n(total points) - n(points in the mean curve)
p_value = 1 - chi2cdf(chisq_res, (size(resid,1)-n_para));