function [Tcorr t cross_prod_binned]= compute_ACF2(x,dt, block_min_ind, varargin) % varargin method options, empty default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute TACF and its covariance matrix using the multi-tau algorithm and
% the block-transformation
%
% ----Required input------------------------------------------------------
% x:  2D array of 8-bit intensity traces with individual traces in rows
% dt:  sampling time of the intensity trace (sec)
% block_min_ind:  index of the minimal block time from block analysis
% ----Optional input------------------------------------------------------
% varargin:  limit for the lag time (sec) 
% ----Output-------------------------------------------------------
% Tcorr:  2D array of TACF curves with individual curves in rows
% t:  column vector of the lag time (sec)
% cross_prod_binned: 3D array of block-averaged intensity products given by 
% the block-transformation with
%   d1-intensity products
%   d2-lag time
%   d3-intensity traces
%-----------------------------------------------------------------
% Ref: 1. T. Wohland et. al., BPJ, 2001
%      2. S. Guo et. al., Anal Chem, 2012
% -----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Apr 02, 2012
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Calculating time correlation functions...');
if ~isempty(varargin)
    limit = varargin{1};
else
    limit = 0.1* size(x,2) ; % default limit for the lag time
end      
   % set semilog scale increment
   increment(1:8)=1;
    for k = 1:30
        increment(k*8+1:(k+1)*8)=2^(k-1);
    end    
    t = cumsum(increment);
    t = t(1:sum(t<limit));
    nt = size(t,2); % number of lags
    x_sqz = squeeze(x) ;
    num_frame = size(x_sqz,2) ; % number of intensity values          
    bin_ct = 2.^(block_min_ind-1);
    Tcorr = zeros(size(x_sqz,1),size(t,2));
    bin_c = max(bin_ct./increment(length(t)),ones(size(bin_ct))) ; % minimum binning is 1    
    n_prod_min = floor(num_frame/increment(length(t)))-t(end)/increment(length(t)) ;  % minimum number of products (at maximum lag time)         
%     cross_prod_binned = zeros(floor(n_prod_min/bin_c),size(t,2), size(x_sqz,1));            
    cross_prod_binned = cell(size(x_sqz,1));
   for m = 1: size(x_sqz,1)
     x1 = x_sqz(m,:) ;
     incre_current = 1;
     x1 = var_expd(x1) ; % check and expand the size of the variable according to the values of the variable
     for j = 1 : nt
        if incre_current ~= increment(j) % bin the intensity trace and change the variable type accordingly to save the memories         
             x1 = binning_1d(x1,2,'native');
             x1 = var_expd(x1) ; % check and expand the size of the variable according to the values of the variable                 
             incre_current = increment(j);
        end
        lag = t(j)/incre_current; % update the lag because the intensity trace is binned 
        mean_del = mean(x1(lag+1:end));
        mean_dir = mean(x1(1:end-lag));                     
        cross_prod = x1(lag+1:end).*x1(1:end-lag) ; 
        Tcorr(m,j) = mean(cross_prod)./mean_del./mean_dir ;
        prod_bsize = increment(length(t))/incre_current ;
        % block-transformation
        cross_prod_block_s = double(binning_1d(cross_prod,prod_bsize* bin_c(m), 'double' ))/(prod_bsize* bin_c(m)) ;
        loc_m_del = double(binning_1d(x1(lag+1:end), prod_bsize* bin_c(m),'double'  ))/(prod_bsize* bin_c(m)) ;
        loc_m_dir = double(binning_1d(x1(1:end-lag), prod_bsize* bin_c(m),'double' ))/(prod_bsize* bin_c(m)) ;        
        % convert <F(t)F(t+tau)> to <dF(t)dF(t+tau)>
        cross_prod_block_s = cross_prod_block_s - mean_del*loc_m_dir ...
            - mean_dir*loc_m_del + mean_del*mean_dir ;
        % size(cross_prod_block_s)
        cross_prod_binned{m}(:,j) = cross_prod_block_s(1:floor(n_prod_min/bin_c(m)))./mean_del./mean_dir ;               
        waitbar((m-1)/size(x_sqz,1)+j/size(t,2)*1/size(x_sqz,1),h)
     end        
   end

    Tcorr = Tcorr - 1;
    close(h)
    t = t'*dt;      
end

function x1 = var_expd(x1) % to save the RAM usage, check and expand the size of the variable according to the values of the variable
    if max(x1(:))>8 && max(x1(:))<=2^7 && ~isa(x1, 'uint16')
        x1 = uint16(x1); 
    elseif max(x1(:))>2^7 && max(x1(:))<= 2^15 && ~isa(x1, 'uint32')
        x1 = uint32(x1); 
    elseif max(x1(:))>2^15 && ~isa(x1, 'uint64')
        x1 = uint64(x1);
    else
    end            
end
