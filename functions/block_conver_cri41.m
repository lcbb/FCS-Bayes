function block_min_ind = block_conver_cri41(ste_it, err)
% This function determines the minimal block-time (fixed point) in blocking 
% curves. It seeks the first three consecutive points in the blocking curve 
% that satisfy the following convergence criteria: 
% Criterion 1. The differences between nearest-neighbor, or consecutive, blocking points in the series must be smaller than the average of their error bars.
% Criterion 2. No three consecutive points follow the series, which are (i) monotonically increasing and (ii) violate criterion 1.
% Criterion 3. The last point in the series does not belong to the ultimate or penultimate point in the blocking curve. 
%-------Input----------------------------------------------------
% ste_it: 2D array of blocking curves with individual curves in rows
% err: error bars associated with blocking curves (standard
% errors).
%-------Ouput----------------------------------------------------
% block_min_ind: row vector of the index of the fixed points for blocking curves
%-----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Oct 31, 2013
%-----------------------------------------------------------------
block_min_ind = zeros(1, size(ste_it,1)) ;
    cr_win = 2 ;    % window size of 3 points (2 pair-wise comparisons)
    for cur = 1:size(ste_it,1)
        ste_it_s = squeeze(ste_it(cur,:));    % take a single curve
        err_s = squeeze(err(cur,:)); 
        ste_diff = diff(ste_it_s) ;
        cr1 = abs(ste_diff)<(err_s(1:end-1)+ err_s(2:end)) ; % pairs of points with their error bars overlapping
        cr1 = tsmovavg((cr1), 's', cr_win)'>=1 ; % 3 points with their error bars overlapping       
        diff_pos = ste_diff>0 ; % pairs with poisitive difference 
        cr2 = tsmovavg(diff_pos, 's', cr_win)'<1 ; % indices of 3 points with at least one of differences being negative                              
        last_0 = find((~cr2&~cr1), 1, 'last') ; % index of the last point violating Criterion 2                              
        cr12 = cr1 ; 
        cr12(1:last_0) = 0 ; % indices of two-pairs that pass Criterion 1 & 2
        cr12(end-1:end)=0 ; % Criterion 3
        ind = find(cr12==1, 1, 'first') ;
        if ~isempty(ind)
            block_min_ind(cur) = ind ; 
        else
        end
    end