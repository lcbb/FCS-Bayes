function x_binned = binning_3d(x,bin_size,varargin)
if ~isempty(varargin)
    opt = varargin{1} ;    
else
    opt = 'double' ;
end

num_new  = floor(size(x,3)/bin_size);
switch opt
    case 'native'        
        x_binned  = zeros(size(x,1), size(x,2), num_new,class(x));
    case 'double'
        x_binned  = zeros(size(x,1), size(x,2), num_new);
end

for j = 1 : num_new
    x_binned(:,:,j) = sum(x(:,:,(j-1)*bin_size+1:j*bin_size),3,opt);
end
end