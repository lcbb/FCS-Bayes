function [block_min_ind block_cur block_err t_block]= blocking_im(im, dt, block_lag, limit, plot_opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function runs the blocking analysis on the movie to
% determine the minimum block-time needed to obtain indepedent samples for
% proper estimation of covariance matrices of resulting TACFs.
% - perform blocking analysis on the input traces at the given lag index
% - plot blocking curves
%
% ----Required input------------------------------------------------------
% x:  3D array of the movie (pixel, pixel, time)
% dt: frame interval of the movie (sec)
% t:  column vector of the lag time
% ----Optional input------------------------------------------------------
% block_lag: the lag index where blocking is performed. Default is 1.
% limit: limit for the lag time (sec) 
% plot_opt: option for plotting the blocking results (1) or not (0)
% ----Output-------------------------------------------------------
% block_min_ind: 2D matrix of the index of the fixed points for blocking curves (pixel,
% pixel) 
% block_cur: 3D array of blocking curves of intensity traces (pixel,
% pixel, block time).
% block_error:  error bars associated with blocking curves (standard
% errors).
% t_block:  vector of the block-time
%
%
%-----------------------------------------------------------------
% Ref: 1. Guo, S. et al., Anal Chem 2012
%      2. Flyvbjerg, H. et al., J. Chem. Phys. 1989
% -----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Oct 31, 2013
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default options
if nargin < 5, plot_opt = 0 ; end % default limit for the lag time (10%)
if nargin < 4, limit = 0.1* size(x,2) ; end % default limit for the lag time (10%)
if nargin < 3, block_lag = 1 ;  end % default 

block_min_ind = zeros(size(im,1), size(im,2)) ; % matrix of the index of the fixed points 
num_frame = size(im,3);
% set fine semilog scale increment
increment(1:16)=1;
    for k = 1:15
        increment(k*16+1:(k+1)*16)=2^(k-1);
    end
n_block_t = floor(log2(num_frame/increment(block_lag)/8)) ; % number of the block transformations performed
t_block = 2.^[0:n_block_t]*increment(block_lag)*dt ; % vector of the block-time
block_cur = zeros(size(im,1), size(im,2), n_block_t+1) ; % blocking curve matrix
block_err = zeros(size(im,1), size(im,2), n_block_t+1) ; % blocking error matrix
t = cumsum(increment); 
t = t(1:sum(t<limit)); % lag time vector for the TACF

h = waitbar(0,'Running blocking analysis...'); 

for m=1:size(im,1)
    x1 = permute(im(m,:,:),[2 3 1]) ;                  
    % select the binning function based on the dimension of the image
    if size(x1,1)== 1          
        binning = @binning_1d ;
    else
        binning = @binning_2d ;
    end   
    
    var_it = zeros(size(x1,1),size(t,2),n_block_t+1);  % variance
    ste_it = zeros(size(x1,1),size(t,2),n_block_t+1);  % standard errors
    n_prod = zeros(size(x1,1),size(t,2),n_block_t+1);             
    incre_current = 1; % currect increment
    
    for j = 1 : block_lag
        % bin the intensity trace using the multi-tau correlator scheme
        if incre_current ~= increment(j)           
             x1 = binning(x1,2,'native');                         
             incre_current = increment(j);
        else
        end
        
        if j == block_lag ; 
            lag = t(j)/incre_current; % update the lag because the intensity trace is binned     
            mean_del = mean(x1(:,lag+1:end),2);
            mean_dir = mean(x1(:,1:end-lag),2);     
            cross_prod = bsxfun(@minus,x1(:,lag+1:end),mean_del).*bsxfun(@minus,x1(:,1:end-lag),mean_dir) ; % dF(t)dF(t+tau)                                                              
            var_it(:,j,1) = var(cross_prod,0,2)/size(cross_prod,2)./mean_del.^2./mean_dir.^2  ;
            n_prod(:,j,1) = size(cross_prod,2) ;          

            for k = 1:n_block_t
                cross_prod = binning_2d(cross_prod,2)/2 ;
                var_it(:,j,k+1) = var(cross_prod,0,2)/size(cross_prod,2)./mean_del.^2./mean_dir.^2 ;  
                n_prod(:,j,k+1) = size(cross_prod,2) ;
            end                     
        end 
    end    
    ste_it(:,block_lag,:) = sqrt(var_it(:,block_lag,:));   % standard errors
    err = ste_it.*1./sqrt(2*(n_prod-1)) ; % errors of standard errors    
    
    % save blocking curves and error
    block_cur(m,:,:) = permute(ste_it(:,block_lag,:), [2 1 3]) ;
    block_err(m,:,:) = permute(err(:,block_lag,:),[2 1 3]) ;    
    ste_it_sb = squeeze(ste_it(:,block_lag,:));    
    err_sb = squeeze(err(:,block_lag,:));    
    block_min_ind(m,:) = block_conver_cri41(ste_it_sb, err_sb) ;  % determine the fixed point using blocking criteria             
    waitbar(m/size(im,1),h)
 % Plot blocking curves 
  if plot_opt ==1
    figure(84)
    plot_blocking(ste_it_sb, err_sb, t_block, block_min_ind(m,:)) 
  else
  end
end
close(h)
end
    