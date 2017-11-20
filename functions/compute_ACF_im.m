function [Tcorr t cross_prod_binned]= compute_ACF_im(x,dt, block_min_ind, varargin) % varargin method options, empty default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute TACFs and their covariance matrices for a movie using the multi-tau algorithm and
% the block-transformation
%
% ----Required input------------------------------------------------------
% x:   3D array of the moive
% dt: frame interval of the movie (sec)
% block_min_ind:  index of the minimal block time from the blocking analysis
% ----Optional input------------------------------------------------------
% varargin:  limit for the lag time (sec) 
% ----Output-------------------------------------------------------
% Tcorr:  3D array of TACF curves (pixel, pixel, lag time).
% t:  column vector of the lag time (sec)
% cross_prod_binned: 2D cell array with each cell recording the matrix of block-averaged intensity products for each pixel.
% The matrix has intensity products in rows and lag time in columns
%-----------------------------------------------------------------
% Ref: 1. T. Wohland et. al., BPJ, 2001
%      2. S. Guo et. al., Anal Chem, 2012
% -----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Oct. 31, 2013
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Calculating time correlation functions...');
num_frame = size(x,3) ;
if ~isempty(varargin)
    limit = varargin{1};
else
    limit = 0.1* size(x,3) ; % default limit for the lag time
end      
   % set fine semilog scale increment
   increment(1:16)=1;
%     for k = 1:15
    for k = 1:9
        increment(k*16+1:(k+1)*16)=2^(k-1);
    end    
    t = cumsum(increment);
    t = t(1:sum(t<limit));        
    
    nt = size(t,2); % number of lags    
    bin_ct = 2.^(block_min_ind-1) ; % binning size given by the minimal block time
    % select the lager binning size of ones given by blocking and multi-tau
    % at maximum tau
    bin_c = max(bin_ct/increment(length(t)),ones(size(bin_ct))) ; % compare the binning size of blocking and the largest channel width
    n_prod_min = floor(num_frame/increment(length(t)))-t(end)/increment(length(t)) ;  % minimum number of products (at maximum lag time) 
    n_prod_binned = floor(n_prod_min./bin_c) ;  % the final number of products after binning    
    Tcorr = zeros(size(x,1),size(x,2),size(t,2)); % TACF
    cross_prod_binned = cell(size(x,1),size(x,2));    
    cov_it = zeros(size(t,2), size(t,2), size(x,1),size(x,2)) ; % covariance matrix for each TACF
    % pre-allocate cross-product cell
    for k = 1:size(x,1)
        for l = 1:size(x,2)                          
            cross_prod_binned{k,l} = zeros(n_prod_binned(k,l), size(t,2)) ;
        end
    end
    
    incre_current = 1;
    for j = 1 : nt 
        % bin the intensity trace using the scheme of the multi-tau
        % correlator
        if incre_current ~= increment(j)
            x = binning_3d(x,2,'native');
            incre_current = increment(j);
        end
        lag = t(j)/incre_current;  % update the lag because the intensity trace is binned 
        mean_del = mean(x(:,:,lag+1:end),3);
        mean_dir = mean(x(:,:,1:end-lag),3);     
        cross_prod = bsxfun(@minus,x(:,:,lag+1:end),mean_del).*bsxfun(@minus,x(:,:,1:end-lag),mean_dir) ; 
%         cross_prod = x(:,:,lag+1:end).*x(:,:,1:end-lag); % F(t)F(t+tau) 
        Tcorr(:,:,j) = mean(cross_prod,3)./mean_del./mean_dir ;
        prod_bsize = increment(length(t))/incre_current ;  
        
        for k = 1:size(x,1)
            for l = 1:size(x,2)                                
                % block-transformation
                cross_prod_block_s = binning_1d(squeeze(cross_prod(k,l,:)), prod_bsize*bin_c(k,l))/(prod_bsize*bin_c(k,l)) ;                
                cross_prod_block_s = cross_prod_block_s(1:n_prod_binned(k,l))/(mean_del(k,l).*mean_dir(k,l)) ;                 
                cross_prod_binned{k,l}(:,j) = cross_prod_block_s ;
            end
        end
        waitbar(j/size(t,2),h) 
    end   
    
    for k = 1:size(x,1)
        for l = 1:size(x,2)            
            cov_it(:,:,k,l) = cov(cross_prod_binned{k,l})./n_prod_binned(k,l)  ;
        end
    end                   
    close(h)
    t = t'*dt;
end

