function [ste_it err t_block]= blocking_analysis_CCF(xa, xb, dt, block_lag, limit)
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
% limit: limit for the lag time (sec) 
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
% -----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Apr 06, 2012
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default options
    if nargin < 4, limit = 0.1* size(xa,2) ; end % default limit for the lag time (10%)
    if nargin < 3, block_lag = 1 ;  end % default 
    
    num_frame = size(xa,2);  % number of intensity values
    h = waitbar(0,'Running blocking analysis...');       
    
    increment(1:8)=1;
    for k = 1:30
        increment(k*8+1:(k+1)*8)=2^(k-1);
    end
    
    % set semilog scale increment
    t = cumsum(increment);
    t = t(1:sum(t<limit));
    nt = size(t,2); % number of lags    
    n_block_t = floor(log2((num_frame-limit)/increment(block_lag)/16)) ; % number of blocking transformations performed                
    var_it = zeros(size(xa,1),n_block_t+1);    % variance
    n_prod = zeros(size(xa,1),n_block_t+1);    % number of blocked samples
                  
     for m = 1: size(xa,1)
        xa1 = xa(m,:) ;
        xb1 = xb(m,:) ;
         incre_current = 1;
         for j = 1 : block_lag
             if incre_current ~= increment(j)
                 xa1 = binning_1d(xa1,2,'native');
                 xb1 = binning_1d(xb1,2, 'native');
                 if min(max(xa1), max(xb1))>8 && min(max(xa1), max(xb1))<=2^7 ...
                     && ~isa(xa1, 'uint16')
                     xa1 = uint16(xa1); xb1 = uint16(xb1);
                 elseif min(max(xa1), max(xb1))>2^7 && min(max(xa1), max(xb1))<= 2^15 ...
                         && ~isa(xa1, 'uint32')
                     xa1 = uint32(xa1); xb1 = uint32(xb1);
                 elseif min(max(xa1), max(xb1))>2^15 && ~isa(xa1, 'uint64')
                     xa1 = uint64(xa1); xb1 = uint64(xb1);
                 end
                 incre_current = increment(j);
            end           
            
            if j == block_lag                
                lag = t(j)/incre_current; % update the lag because the intensity trace is binned             
                mean_del = mean(xa1(lag+1:end));
                mean_dir = mean(xb1(1:end-lag));            
                cross_prod = xa1(lag+1:end).*xb1(1:end-lag) ;  
                % convert <F(t)F(t+tau)> to <dF(t)dF(t+tau)>
                cross_prod = double(cross_prod) - mean_del*double(xb1(1:end-lag)) ...
                    - mean_dir*double(xa1(lag+1:end)) + mean_del*mean_dir ;
                var_it(m,1) = var(cross_prod,0,2)/size(cross_prod,2)./mean_del.^2./mean_dir.^2  ;
                n_prod(m,1) = size(cross_prod,2) ;                
                for k = 1:n_block_t
                    cross_prod = binning_1d(cross_prod,2)/2 ;
                    var_it(m,k+1) = var(cross_prod,0,2)/size(cross_prod,2)./mean_del.^2./mean_dir.^2  ;
                    n_prod(m,k+1) = size(cross_prod,2) ;
                    waitbar((m-1)/size(xa1,1)+k/n_block_t*1/size(xa1,1),h)
                end                
                t_block = 2.^[0:n_block_t]*increment(j)*dt ; % vector of the block-time
            end            
         end
     end
         t = t'*dt;
         ste_it = sqrt(var_it) ;    % standard errors
         err = ste_it.*1./sqrt(2*(n_prod-1)) ;   % errors of standard errors
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