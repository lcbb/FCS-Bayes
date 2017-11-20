% Same as the function "compute_ACF" but calculate the cross-correlation
% funtion. Refer to "compute_ACF" for details.
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function [Tcorr t cross_prod_binned]= compute_CCF(xa,xb,dt,bin_c, varargin)

h = waitbar(0,'Calculating time correlation functions...');
if ~isempty(varargin)
    limit = varargin{1};
else
    limit = 0.1* size(xa,2) ; % default limit for the lag time
end          
     % set semilog scale increment
   increment(1:8)=1;
    for k = 1:20
        increment(k*8+1:(k+1)*8)=2^(k-1);
    end    
    t = cumsum(increment);
    t = t(1:sum(t<limit));
    nt = size(t,2);
    num_frame = size(xa,2) ;
    Tcorr = zeros(size(xa,1),size(t,2));
    n_prod_min = floor(num_frame/increment(length(t)))-t(end)/increment(length(t)) ;  % minimum number of products      
    cross_prod_binned = zeros(floor(n_prod_min/bin_c),size(t,2), size(xa,1));            
     for m = 1: size(xa,1)       
        xa1 = xa(m,:) ;
        xb1 = xb(m,:) ;
         incre_current = 1;
        for j = 1 : nt
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
            lag = t(j)/incre_current; % update the lag because the intensity trace is binned 
            xa_par = zeros(1,size(xa1,2)-lag,class(xa1));
            xb_par = zeros(1,size(xb1,2)-lag,class(xb1));            
            xa_par = xa1(lag+1:end) ;
            mean_del = mean(xa_par);
            xb_par = xb1(1:end-lag) ;
            mean_dir = mean(xb_par);                     
            cross_prod = xa_par.*xb_par ; 
            Tcorr(m,j) = mean(cross_prod)./mean_del./mean_dir ;
            prod_bsize = increment(end)/incre_current ;               
            cross_prod_block_s = double(binning_1d(cross_prod,prod_bsize* bin_c, 'double' ))/(prod_bsize* bin_c) ;
            loc_m_del = double(binning_1d(xa1(lag+1:end), prod_bsize* bin_c,'double'  ))/(prod_bsize* bin_c) ;
            loc_m_dir = double(binning_1d(xb1(1:end-lag), prod_bsize* bin_c,'double' ))/(prod_bsize* bin_c) ;
            % convert <F(t)F(t+tau)> to <dF(t)dF(t+tau)>
            cross_prod_block_s = cross_prod_block_s - mean_del*loc_m_dir ...
                - mean_dir*loc_m_del + mean_del*mean_dir ;
            size(cross_prod_block_s)
            cross_prod_binned(:,j,m) = cross_prod_block_s(1:floor(n_prod_min/bin_c))./mean_del./mean_dir ;            
            waitbar((m-1)/size(xa,1)+j/size(t,2)*1/size(xa,1),h)
        end            
     end

        Tcorr = Tcorr - 1;
        close(h)
        t = t'*dt;
end

