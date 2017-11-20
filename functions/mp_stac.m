% This function calculate quartiles of model probabiliies
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function [mp_m err_l err_h] = mp_stac(mp)
n_condi = numel(mp) ;
n_model = size(mp{1},2) ;
mp_m = zeros(n_model,n_condi); mp_lq = mp_m; mp_hq = mp_m ;
for j = 1:n_condi
    for m = 1:n_model        
        mp_m(m,j) = median(mp{j}(:,m));
        mp_lq(m,j) = prctile(mp{j}(:,m),25);
        mp_hq(m,j) = prctile(mp{j}(:,m),75);
    end
end
err_l = (mp_m - mp_lq) ;
err_l(err_l<0)= 0 ;
err_h = (mp_hq - mp_m) ;
err_h(err_h<0)= 0 ;
end