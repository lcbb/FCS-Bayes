% This funtion bins the each row of the input 2D matrix with "bin_size"
% "opt" specifies the variable type
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function x_binned = binning_2d(x,bin_size,varargin)

if ~isempty(varargin)
    opt = varargin{1} ;    
else
    opt = 'double' ;
end

num_new  = floor(size(x,2)/bin_size);
switch opt
    case 'native'        
        x_binned  = zeros(size(x,1),num_new,class(x));
    case 'double'
        x_binned  = zeros(size(x,1),num_new);
end

for j = 1 : num_new
    x_binned(:,j) = sum(x(:,(j-1)*bin_size+1:j*bin_size),2,opt);
end
end