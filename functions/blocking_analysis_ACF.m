function [ste_it err t_block n_prod]= blocking_analysis_ACF(x, dt, block_lag, limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function runs the blocking analysis on photon count traces to
% determine the minimum block-time needed to obtain indepedent samples for
% proper estimation of covariance matrices of resulting TACFs.
% - perform blocking analysis on the input traces at the given lag index
% - plot blocking curves
%
% ----Required input------------------------------------------------------
% x:  2D array of 8-bit intensity traces with individual traces in rows
% dt:  sampling time of the intensity trace
% t:  column vector of the lag time
% ----Optional input------------------------------------------------------
% block_lag: the lag index where blocking is performed. Default is 1.
% limit: limit for the lag time in the TACF (sec) 
% ----Output-------------------------------------------------------
% ste_it:  2D array of standard errors estimated from samples with
% different block-time (blocking curves). Each row is the result of one trace.
% t_block:  vector of the block-time
% n_prod:  2D array of numbers of underlying samples for ste_it
%
%
%-----------------------------------------------------------------
% Ref: 1. Guo, S. et al., Anal Chem 2012
%      2. Flyvbjerg, H. et al., J. Chem. Phys. 1989
%      3. T. Wohland et. al., BPJ, 2001
% -----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Apr 06, 2012
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default options
    if nargin < 4, limit = 0.1* size(x,2) ; end % Use 10% of the measurement time as the default lag time limit if not specified
    if nargin < 3, block_lag = 1 ;  end % Run the blocking analysis on intensity products with minimal lag time if not specified 
    
    num_frame = size(x,2);  % number of intensity values in the trace
    h = waitbar(0,'Running blocking analysis...');  % display the wait bar     
    
    % set semilog scale increment for TACF lag time
    increment(1:8)=1; % 8 channels in the first group
    for k = 1:30 % total 31 groups
        increment(k*8+1:(k+1)*8)=2^(k-1);
    end    
    
    t = cumsum(increment); % vector of lags
    t = t(1:sum(t<limit)); % make maximal lag < "limit"
    nt = size(t,2); % number of lags    
    n_block_t = floor(log2((num_frame-limit)/increment(block_lag)/8)) ; % calculate number of blocking transformations to perform to make the number of blocks >8 at the maximal block time                 
    var_it = zeros(size(x,1),n_block_t+1);    % pre-allocate the variable for the variance (Skk in ref.1) at different block-times
    n_prod = zeros(size(x,1),n_block_t+1);    % pre-allocate the variable for the number of blocks (M' in ref.1) at different block-times
                  
     for m = 1: size(x,1) % loop through each trace
         x1 = x(m,:) ;
         incre_current = 1; % currect increment
         for j = 1 : block_lag % loop through each lag
            if incre_current ~= increment(j) % bin the intensity trace according to the increment (channel width), and change the variable type accordingly to save the memory usage         
                 x1 = binning_1d(x1,2,'native'); % bin the intensity trace with binning size 2
                 if max(x1(:))>8 && max(x1(:))<=2^7 && ~isa(x1, 'uint16') 
                     x1 = uint16(x1); % make "x1" 16 bit integer
                 elseif max(x1(:))>2^7 && max(x1(:))<= 2^15 && ~isa(x1, 'uint32')
                     x1 = uint32(x1);  % make "x1" 32 bit integer 
                 elseif max(x1(:))>2^15 && ~isa(x1, 'uint64') 
                     x1 = uint64(x1); % make "x1" 32 bit integer
                 end            
                 incre_current = increment(j); % update the current increment (channel width)
            end           
            
            if j == block_lag  % if lag == block_lag, run blocking analysis             
                lag = t(j)/incre_current; % update the lag according the current increment              
                mean_del = mean(x1(lag+1:end)); % "Mdel" in Ref. 3
                mean_dir = mean(x1(1:end-lag)); % "Mdir" in Ref. 3           
                cross_prod = x1(lag+1:end).*x1(1:end-lag) ; % F(t)F(t+tau)                
                cross_prod = double(cross_prod) - mean_del*double(x1(1:end-lag)) ...
                    - mean_dir*double(x1(lag+1:end)) + mean_del*mean_dir ; % convert <F(t)F(t+tau)> to <dF(t)dF(t+tau)>
                var_it(m,1) = var(cross_prod,0,2)/size(cross_prod,2)./mean_del.^2./mean_dir.^2  ; % Eq. 13 in Ref. 1
                n_prod(m,1) = size(cross_prod,2) ;                
                for k = 1:n_block_t
                    cross_prod = binning_1d(cross_prod,2)/2 ; % blocking transformation (Eq. 20 in Ref. 2, Eq. 12 in Ref. 1)
                    var_it(m,k+1) = var(cross_prod,0,2)/size(cross_prod,2)./mean_del.^2./mean_dir.^2  ; % Var(<dF(t)dF(t+tau)>/<F>^2) = Var(dF(t)dF(t+tau))/N/<F>^2
                    n_prod(m,k+1) = size(cross_prod,2) ;
                    waitbar((m-1)/size(x,1)+k/n_block_t*1/size(x,1),h) % update the wait bar
                end                
                t_block = 2.^[0:n_block_t]*increment(j)*dt ; % vector of the block-time
            end            
         end
     end
         t = t'*dt; % lag time =  lags * sampling time
         ste_it = sqrt(var_it) ;    % standard errors
         err = ste_it.*1./sqrt(2*(n_prod-1)) ;   % errors of standard errors (Eq. 28 in Ref. 2)
    close(h)
    
%% Plot blocking curves
% figure(8)
% for cur = 1:size(ste_it,1)
%     errorbar(t_block, squeeze(ste_it(cur,:)),squeeze(err(cur,:)),'linewidth',2, 'markersize',8);
%     hold all
% end
% hold off
% set(gca,'xscale','log');
% set(gca, 'xtick',10.^[-6:6])
% axis tight
% xlim = get(gca, 'xlim') ; ylim = get(gca, 'ylim') ;
% % set(gca, 'ylim', [0 ylim(2)])
% title(['\tau = ' num2str(t(block_lag)) ' s'], 'fontsize', 20) ;
% ylabel('\sigma','FontSize',10)
% xlabel('block time (s)','FontSize',10)
% format_fig2(5)

end